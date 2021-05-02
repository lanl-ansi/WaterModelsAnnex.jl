import Random

function get_controllable_variable_indices(wm::WM.AbstractWaterModel, n::Int)
    var_ids = Array{WM._VariableIndex, 1}([])

    for comp_type in [:pump, :valve, :regulator]
        comp_ids = sort(collect(WM.ids(wm, n, comp_type)))
        var_sym = Symbol("z_" * String(comp_type))
        append!(var_ids, [WM._VariableIndex(n, comp_type, var_sym, i) for i in comp_ids])
    end

    return var_ids
end


function compare_variable_indices(var_1::WM._VariableIndex, var_2::WM._VariableIndex)
    network_index_match = var_1.network_index == var_2.network_index
    component_type_match = var_1.component_type == var_2.component_type
    variable_symbol_match = var_1.variable_symbol == var_2.variable_symbol
    component_index_match = var_1.component_index == var_2.component_index
    return all([network_index_match, component_type_match, variable_symbol_match, component_index_match])
end


mutable struct ControlSetting
    network_id::Int
    variables_indices::Array{WM._VariableIndex, 1}
    values::Array{Float64, 1}
end


function create_pump_settings_at_time(wm::WM.AbstractWaterModel, n::Int)
    control_settings = Array{Array{ControlSetting}}([])

    # Create variable indices for all controllable components at time `n`.
    for (pump_group_id, pump_group) in WM.ref(wm, n, :pump_group)
        pump_ids = pump_group["pump_indices"]
        pump_combinations = vcat([zeros(length(pump_ids))], [[j >= i ? 1.0 :
            0.0 for i = 1:length(pump_ids)] for j = 1:length(pump_ids)])
        pump_vars = WM._VariableIndex.(n, :pump, :z_pump, pump_ids)
        push!(control_settings, ControlSetting.(n, Ref(pump_vars), pump_combinations))
    end

    pump_ids_group = vcat([collect(x["pump_indices"])
        for (i, x) in WM.ref(wm, n, :pump_group)]...)

    pump_ids = filter(x -> !(x in pump_ids_group), WM.ids(wm, n, :pump))

    # if length(pump_ids) > 0
    #     z = WM._VariableIndex.(n, :pump, :z_pump, pump_ids)
    #     #vals = collect(Iterators.product([1.0, 0.0] for k in 1:length(z)))
    #     vals = vcat(collect(Iterators.product([[1.0, 0.0] for k in 1:length(z)])...))
    #     push!(control_settings, ControlSetting.(n, Ref(z), vals))
    # end

    return control_settings
end


function create_regulator_settings_at_time(wm::WM.AbstractWaterModel, n::Int)
    regulator_ids = WM.ids(wm, n, :regulator)
    control_settings = Array{Array{ControlSetting}}([])

    if length(regulator_ids) > 0
        z = WM._VariableIndex.(n, :regulator, :z_regulator, regulator_ids)
        vals = collect(Iterators.product([1.0, 0.0] for k in 1:length(z)))
        push!(control_settings, ControlSetting.(n, Ref(z), vals))
    end

    return control_settings
end


function create_valve_settings_at_time(wm::WM.AbstractWaterModel, n::Int)
    valve_ids = WM.ids(wm, n, :valve)
    control_settings = Array{Array{ControlSetting}}([])

    if length(valve_ids) > 0
        z = WM._VariableIndex.(n, :valve, :z_valve, valve_ids)
        vals = collect(Iterators.product([1.0, 0.0] for k in 1:length(z)))
        push!(control_settings, ControlSetting.(n, Ref(z), vals))
    end

    return control_settings
end


function cartesian_product(settings::Array{Array{ControlSetting}})
    return collect(Iterators.product(settings...))
end


function create_control_settings_at_time(wm::WM.AbstractWaterModel, n::Int)
    pump_control_settings = create_pump_settings_at_time(wm, n)
    # regulator_control_settings = create_regulator_settings_at_time(wm, n)
    # valve_control_settings = create_valve_settings_at_time(wm, n)
    # control_settings_partitioned = vcat(pump_control_settings,
    #     regulator_control_settings, valve_control_settings)

    control_settings = cartesian_product(pump_control_settings)
    

    

    # z_valve = WM._VariableIndex.(n, :valve, :z_valve, WM.ids(wm, n, :valve))
    # append!(control_settings, [ControlSetting.(n, [[v] for v in z_valve], Ref([1.0]))])
    # append!(control_settings, [ControlSetting.(n, [[v] for v in z_valve], Ref([0.0]))])

    # Generate all possible controllable status combinations at time `n`.
    # var_ids = vcat(z_pump, z_regulator, z_valve)
    # combs = Iterators.product([[0.0, 1.0] for k in 1:length(var_ids)]...)
    # append!(control_settings, [ControlSetting(n, var_ids, collect(c)) for c in combs])
    return pump_control_settings
    #return vcat([ControlSetting(n, var_ids, collect(c)) for c in combinations]...)
end


function create_all_control_settings(wm::WM.AbstractWaterModel)
    network_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]
    return vcat([create_control_settings_at_time(wm, n) for n in network_ids]...)
end


function filter_control_settings_pump_group!(control_settings::Array{ControlSetting, 1})
    println("hi")
end

    #return vcat([create_control_settings_at_time(Ref(wm), n) for n in network_ids]...)

    # schedules = Dict{Int, Tuple}()

    # for n in sort(collect(WM.nw_ids(wm)))[1:end-1]
    #     # Get all controllable components at time step `n`.
    #     var_ids = get_controllable_variable_indices(wm, n)

    #     # Create and store all possible component schedules at time step `n`.
    #     combinations = collect(Iterators.product([[0, 1] for k in 1:length(var_ids)]...))
    #     schedules[n] = (var_ids, combinations)

    #     for (k, pump_group) in WM.ref(wm, n, :pump_group)
    #         pump_ids = sort(collect(pump_group["pump_indices"]))
    #         var_seq = [WM._VariableIndex(n, :pump, :z_pump, i) for i in pump_ids]
    #         var_indices = [findfirst(x -> compare_variable_indices(x, var), schedules[n][1]) for var in var_seq]
    #         new = filter(x -> sort(collect(x[var_indices]); rev = true) == collect(x[var_indices]), schedules[n][2])
    #         schedules[n] = (var_ids, collect(new))
    #     end
    # end

    # return schedules



function simulate_schedules(data::Dict{String,<:Any}, schedules::Tuple, optimizer::Any)
    schedule_result = Array{Tuple}([])

    for schedule in schedules[2]
        for (i, var) in enumerate(schedules[1])
            comp_type = string(var.component_type)
            comp_index = string(var.component_index)
            comp = data[comp_type][comp_index]
            comp["status"] = WM.STATUS(schedule[i])
        end

        WM.fix_all_indicators!(data)
        result = WM.solve_wf(data, CDWaterModel,
            optimizer; relax_integrality = true)

        if !feasible_simulation_result(result)
            result["primal_status"] = INFEASIBLE_POINT
            result["objective"] = Inf
        else
            result["primal_status"] = FEASIBLE_POINT
            result["objective"] = 0.0
        end

        solution_has_pumps = haskey(result["solution"], "pump")

        if solution_has_pumps && length(result["solution"]["pump"]) > 0
            cost = sum(pump["c"] for (a, pump) in result["solution"]["pump"])
        else
            cost = 0.0 # If there are no pumps, the cost must be zero.
        end

        tank_sol = Dict{String, Any}(i => t["q"] for (i, t) in result["solution"]["tank"])
        sol = Dict{String, Any}("cost" => cost, "tank" => tank_sol)
        push!(schedule_result, (result["primal_status"], sol))
    end

    return schedule_result
end


function _instantiate_cq_model(data::Dict{String, <:Any})
    return WM.instantiate_model(data, CQWaterModel, WM.build_wf)
end


function calc_control_settings(wm::WM.AbstractWaterModel)
end


function calc_possible_schedules(network::Dict{String, <:Any}, mip_optimizer::Any, nlp_optimizer::Any)
    control_settings = create_all_control_settings(network)

    network_median = deepcopy(network)
    _set_median_tank_time_series!(network_median)

    # Create the multinetwork version of the network.
    network_mn = WM.make_multinetwork(network_median)

    # Solve a relaxation of the OWF to get an estimate of tank volume time series.
    ext = Dict(:pipe_breakpoints => 25, :pump_breakpoints => 25)
    wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf; ext = ext)
    result_mip = WM.optimize_model!(wm; optimizer = nlp_optimizer, relax_integrality = true)

    # Get all component schedules across all possible time steps.
    schedules = create_all_schedules(wm)
    return schedules

    # Create copies of the network to compute everything in parallel.
    network_tmp = [deepcopy(network_median) for i in 1:Threads.nthreads()]
    wms = [_instantiate_cq_model(network_tmp[i]) for i in 1:Threads.nthreads()]
    return wms

    # time = @elapsed Threads.@threads for n in sort(collect(keys(schedules)))
    #     # Simulate schedules and filter the solutions.
    #     WM._IM.load_timepoint!(network_tmp[Threads.threadid()], n)
    #     results = simulate_schedules(network_tmp[Threads.threadid()], schedules[n], nlp_optimizer)
    #     ids = filter(k -> results[k][1] === FEASIBLE_POINT, 1:length(results))
    #     schedules[n] = (schedules[n][1], [(schedules[n][2][i], results[i][2]) for i in ids])
    # end

    # # Return schedules with relevant solution properties.
    # return schedules
end


function solve_heuristic_problem(network::Dict{String, <:Any}, schedules, optimizer)
    # Create the multinetwork version of the network.
    network_mn = WM.make_multinetwork(network);
    ref, model = WM.build_ref(network_mn), JuMP.Model()
    nw_ids = sort(collect(keys(ref[:it][WM.wm_it_sym][:nw])))

    # Initialize the objective.
    cost = JuMP.AffExpr(0.0)

    # Initialize the tank head variables
    h = Dict{Int, Any}(n => JuMP.@variable(
        model, [i in keys(ref[:it][WM.wm_it_sym][:nw][n][:tank])],
        base_name = "h[$(n)]") for n in nw_ids[1:end])

    lambda = Dict{Int, Any}(n => JuMP.@variable(
        model, [k in 1:length(schedules[n][2])],
        lower_bound = 0.0, upper_bound = 1.0, start = 0.0,
        base_name = "lambda[$(n)]") for n in nw_ids[1:end-1])

    for n in nw_ids
        time_step = ref[:it][WM.wm_it_sym][:nw][n][:time_step]

        for (i, tank) in ref[:it][WM.wm_it_sym][:nw][n][:tank]
            surface_area = 0.25 * pi * tank["diameter"]^2
            node = ref[:it][WM.wm_it_sym][:nw][n][:node][tank["node"]]

            if n == 1
                head = node["elevation"] + tank["init_level"]
                JuMP.set_lower_bound(h[n][i], head)
                JuMP.set_upper_bound(h[n][i], head)
            end

            if n < nw_ids[end]
                @assert length(lambda[n]) > 0
                
                cost += sum(lambda[n][k] *
                    schedules[n][2][k][2]["cost"]
                    for k in 1:length(lambda[n]))

                JuMP.@constraint(model, sum(lambda[n]) == 1.0)
                dh_sum = JuMP.AffExpr(0.0)

                for k in 1:length(lambda[n])
                    q_tank = schedules[n][2][k][2]["tank"][string(i)]
                    dh_sum += q_tank * lambda[n][k] / surface_area * time_step
                end

                JuMP.@constraint(model, h[n][i] - dh_sum == h[n+1][i])
                JuMP.set_lower_bound(h[n][i], node["head_min"])
                JuMP.set_upper_bound(h[n][i], node["head_max"])
            end

            if n == nw_ids[end]
                head = node["elevation"] + tank["init_level"]
                JuMP.set_lower_bound(h[n][i], head)
                JuMP.set_upper_bound(h[n][i], node["head_max"])
            end
        end
    end

    JuMP.@objective(model, MOI.MIN_SENSE, cost)
    JuMP.set_optimizer(model, optimizer)
    JuMP.optimize!(model)

    if JuMP.primal_status(model) === MOI.FEASIBLE_POINT
        return Dict{Int, Any}(n => JuMP.value.(lambda[n]) for n in nw_ids[1:end-1])
    else
        return nothing
    end
end


function solve_heuristic_master(network::Dict{String, <:Any}, schedules, weights, nlp_optimizer, optimizer)
    # Create the multinetwork version of the network.
    network_mn = WM.make_multinetwork(network)
    result_mn = WM.solve_mn_owf(network_mn, WM.LRDWaterModel,
        optimizer; relax_integrality = true)
    Random.seed!(0)

    ref = WM.build_ref(network_mn)
    nw_ids = sort(collect(keys(ref[:it][WM.wm_it_sym][:nw])))
    delta_star = Dict{Int, Array}(n => zeros(length(schedules[n][1])) for n in nw_ids[1:end-1])

    for n in nw_ids[1:end-1]
        for (i, var_index) in enumerate(schedules[n][1])
            for (k, schedule) in enumerate(schedules[n][2])
                delta_star[n][i] += schedule[1][i] * weights[n][k]
            end
        end
    end

    feasible, num_iterations = false, 0
    previous_sols = Array{Dict{Int, Any}}([])
    n_min = Array{Int64}([])

    while !feasible && num_iterations <= 100
        model = JuMP.Model()
        cost_1 = JuMP.@expression(model, 0.0)
        costs_2 = [JuMP.AffExpr(0.0) for i in 1:length(schedules[1][1])]

        z = Dict{Int, Any}(n => JuMP.@variable(
            model, [k in 1:length(schedules[n][1])],
            binary = true) for n in nw_ids[1:end-1])

        for (k, previous_sol) in enumerate(previous_sols)
            one_vars = Array{JuMP.VariableRef, 1}([])
            zero_vars = Array{JuMP.VariableRef, 1}([])

            for n in nw_ids[1:end-1]
                if n_min[k] !== nothing && n <= n_min[k]
                    one_ids = findall(x -> isapprox(x, 1.0; atol = 1.0e-4), previous_sol[n])
                    zero_ids = findall(x -> isapprox(x, 0.0; atol = 1.0e-4), previous_sol[n])
                    append!(one_vars, [z[n][j] for j in one_ids])
                    append!(zero_vars, [z[n][j] for j in zero_ids])
                end
            end

            JuMP.@constraint(model, sum(zero_vars) - sum(one_vars) >= 1.0 - length(one_vars))
        end

        for (i, var_index) in enumerate(schedules[1][1])
            for n in nw_ids[1:end-1]
                if isapprox(delta_star[n][i], 0.0; atol = 1.0e-4)
                    Delta = 0.0
                elseif isapprox(delta_star[n][i], 1.0; atol = 1.0e-4)
                    Delta = 1.0
                elseif num_iterations == 0
                    Delta = delta_star[n][i]
                else
                    Delta = Random.shuffle([0.0, delta_star[n][i], 1.0])[1]
                end

                cost_1 += (Delta - z[n][i])^2
                costs_2[i] += Delta - z[n][i]
            end
        end

        cost_2 = sum([costs_2[i]^2 for i in 1:length(schedules[1][1])])
        JuMP.@objective(model, MOI.MIN_SENSE, cost_1 + cost_2)
        JuMP.set_optimizer(model, optimizer)
        JuMP.optimize!(model)
        num_iterations += 1

        @assert JuMP.primal_status(model) === FEASIBLE_POINT        

        sol = Dict{Int, Any}(n => [round(JuMP.value(z[n][k])) for
            k in 1:length(schedules[1][1])] for n in nw_ids[1:end-1])

        simulate!(network, schedules, sol, result_mn, nlp_optimizer)
        feasible = feasible_simulation_result(result_mn)

        if !feasible
            result_mn["primal_status"] = INFEASIBLE_POINT
            result_mn["objective"] = Inf
        else
            result_mn["primal_status"] = FEASIBLE_POINT
            result_mn["objective"] = 0.0
        end

        if result_mn["last_nw"] !== nothing
            sol_reduced = Dict{Int, Any}(n => [round(JuMP.value(z[n][k])) for
                k in 1:length(schedules[1][1])] for n in 1:result_mn["last_nw"])

            push!(n_min, result_mn["last_nw"])
            push!(previous_sols, sol_reduced)
        elseif feasible
            delete!(result_mn, "last_nw")
        end
    end

    return result_mn
end
import Random

function get_controllable_variable_indices(wm::WM.AbstractWaterModel, n::Int)
    var_ids = Array{WM._VariableIndex, 1}([])

    for comp_type in [:pump, :valve, :regulator]
        comp_ids = collect(WM.ids(wm, n, comp_type))
        var_sym = Symbol("z_" * String(comp_type))
        append!(var_ids, [WM._VariableIndex(n, comp_type, var_sym, i) for i in comp_ids])
    end

    return var_ids
end


function create_all_schedules(wm::WM.AbstractWaterModel)
    schedules = Dict{Int, Tuple}()

    for n in WM.nw_ids(wm)
        # Get all controllable components at time step `n`.
        var_ids = get_controllable_variable_indices(wm, n)

        # Create and store all possible component schedules at time step `n`.
        combinations = collect(Iterators.product([[0, 1] for k in 1:length(var_ids)]...))
        schedules[n] = (var_ids, combinations)
    end

    return schedules
end


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
        result = WM.solve_wf(data, CDWaterModel, optimizer; relax_integrality = true)

        objective_is_small = result["objective"] <= 1.0e-6
        status_feasible = result["primal_status"] === FEASIBLE_POINT

        if !(objective_is_small && status_feasible)
            result["primal_status"] = INFEASIBLE_POINT
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


function calc_possible_schedules(network::Dict{String, <:Any}, mip_optimizer::Any, nlp_optimizer::Any)
    # Create the multinetwork version of the network.
    network_mn = WM.make_multinetwork(network)

    # Solve a relaxation of the OWF to get an estimate of tank volume time series.
    ext = Dict(:pipe_breakpoints => 10, :pump_breakpoints => 10)
    wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf; ext = ext)
    result_mip = WM.optimize_model!(wm; optimizer = nlp_optimizer, relax_integrality = true)

    # Get all component schedules across all possible time steps.
    schedules = create_all_schedules(wm)

    # Create copies of the network to compute everything in parallel.
    network_tmp = [deepcopy(network) for i in 1:Threads.nthreads()]

    time = @elapsed Threads.@threads for n in sort(collect(keys(schedules)))
        # Simulate schedules and filter the solutions.
        WM._IM.load_timepoint!(network_tmp[Threads.threadid()], n)
        results = simulate_schedules(network_tmp[Threads.threadid()], schedules[n], nlp_optimizer)
        ids = filter(k -> results[k][1] === FEASIBLE_POINT, 1:length(results))
        schedules[n] = (schedules[n][1], [(schedules[n][2][i], results[i][2]) for i in ids])
    end

    # Return schedules with relevant solution properties.
    return schedules
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
        base_name = "h[$(n)]") for n in nw_ids)

    lambda = Dict{Int, Any}(n => JuMP.@variable(
        model, [k in 1:length(schedules[n][2])],
        lower_bound = 0.0, upper_bound = 1.0, start = 1.0,
        base_name = "lambda[$(n)]") for n in nw_ids)

    for n in nw_ids
        JuMP.@constraint(model, sum(lambda[n]) == 1.0)
        time_step = ref[:it][WM.wm_it_sym][:nw][n][:time_step]
        cost += sum(lambda[n][k] * schedules[n][2][k][2]["cost"] for k in 1:length(lambda[n]))

        for (i, tank) in ref[:it][WM.wm_it_sym][:nw][n][:tank]
            surface_area = 0.25 * pi * tank["diameter"]^2
            node = ref[:it][WM.wm_it_sym][:nw][n][:node][tank["node"]]

            if n == 1
                head = node["elevation"] + tank["init_level"]
                JuMP.set_lower_bound(h[n][i], head)
                JuMP.set_upper_bound(h[n][i], head)
            end

            if n < nw_ids[end]
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

    return Dict{Int, Any}(n => JuMP.value.(lambda[n]) for n in nw_ids)
end


function solve_heuristic_master(network::Dict{String, <:Any}, schedules, weights, nlp_optimizer, optimizer)
    # Create the multinetwork version of the network.
    network_mn = WM.make_multinetwork(network)
    result_mn = WM.solve_mn_owf(network_mn, WM.LRDWaterModel, optimizer; relax_integrality = true)
    Random.seed!(0)

    ref = WM.build_ref(network_mn)
    nw_ids = sort(collect(keys(ref[:it][WM.wm_it_sym][:nw])))
    delta_star = Dict{Int, Array}(n => zeros(length(schedules[n][1])) for n in nw_ids)

    for n in nw_ids
        for (i, var_index) in enumerate(schedules[n][1])
            for (k, schedule) in enumerate(schedules[n][2])
                delta_star[n][i] += schedule[1][i] * weights[n][k]
            end
        end
    end

    feasible, num_iterations = false, 0
    previous_sols = Array{Dict{Int, Any}}([])
    n_min = Array{Int64}([])

    while !feasible && num_iterations <= 1000
        model = JuMP.Model()
        cost_1 = JuMP.@expression(model, 0.0)
        costs_2 = [JuMP.AffExpr(0.0) for i in 1:length(schedules[1][1])]

        z = Dict{Int, Any}(n => JuMP.@variable(
            model, [k in 1:length(schedules[n][1])],
            binary = true) for n in nw_ids)

        for (k, previous_sol) in enumerate(previous_sols)
            one_vars = Array{JuMP.VariableRef, 1}([])
            zero_vars = Array{JuMP.VariableRef, 1}([])

            for n in nw_ids
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
            for n in nw_ids
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

        sol = Dict{Int, Any}(n => [round(JuMP.value(z[n][k])) for
            k in 1:length(schedules[1][1])] for n in nw_ids)

        simulate!(network, schedules, sol, result_mn, nlp_optimizer)
        feasible = result_mn["primal_status"] === FEASIBLE_POINT

        if result_mn["last_nw"] !== nothing
            push!(n_min, result_mn["last_nw"])
        elseif feasible
            delete!(result_mn, "last_nw")
        end
    end

    return result_mn
end


# function _update_control_time_series!(data::Dict{String, <:Any}, schedule::Tuple)
#     if !haskey(data["time_series"], "pump")
#         data["time_series"]["pump"] = Dict{String, Any}()
#     end

#     if !haskey(data["time_series"], "valve")
#         data["time_series"]["valve"] = Dict{String, Any}()
#     end

#     nw_ids = sort([parse(Int, x) for x in keys(result["solution"]["nw"])])

#     for (i, pump) in data["pump"]
#         sol = result["solution"]["nw"]

#         if !haskey(data["time_series"]["pump"], i)
#             data["time_series"]["pump"][i] = Dict{String, Any}()
#         end

#         statuses = Array{Int64, 1}([Int(round(sol[string(n)]["pump"][i]["status"])) for n in nw_ids])
#         data["time_series"]["pump"][i]["status"] = Array{Int64, 1}(statuses)
#     end

#     for (i, valve) in data["valve"]
#         sol = result["solution"]["nw"]

#         if !haskey(data["time_series"]["valve"], i)
#             data["time_series"]["valve"][i] = Dict{String, Any}()
#         end

#         statuses = Array{Int64, 1}([Int(round(sol[string(n)]["valve"][i]["status"])) for n in nw_ids])
#         data["time_series"]["valve"][i]["status"] = statuses
#     end
# end
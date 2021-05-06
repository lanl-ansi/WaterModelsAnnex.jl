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


function set_parameters_from_control_setting!(wm::AbstractCQModel, control_setting::ControlSetting)
    for (i, v) in enumerate(control_setting.variable_indices)
        symbol, id = v.variable_symbol, v.component_index
        JuMP.set_value(WM.var(wm, symbol)[id], control_setting.vals[i])
    end
end


function simulate_control_settings(network::Dict{String, <:Any}, control_settings::Array{ControlSetting, 1}, optimizer)
    networks = [deepcopy(network) for i in 1:Threads.nthreads()]
    wms = [_instantiate_cq_model(networks[i], optimizer) for i in 1:Threads.nthreads()]
    results = [SimulationResult(false, Dict{Int, Float64}(), 0.0) for n in 1:length(control_settings)]
    
    Threads.@threads for k in 1:length(control_settings)
        results[k] = simulate_control_setting(wms[Threads.threadid()], control_settings[k])
    end

    return results
end


function calc_simulation_tank_flows(wm::AbstractCQModel)::Dict{Int, Float64}
    return Dict(i => JuMP.value(WM.var(wm, :q_tank, i)) for i in WM.ids(wm, :tank))
end


function calc_pump_cost(wm::AbstractCQModel, pump_index::Int64)::Float64
    pump = WM.ref(wm, :pump, pump_index)
    q = max(0.0, JuMP.value(WM.var(wm, :q_pump, pump_index)))
    z = JuMP.value(wm.model[:z_pump][pump_index])
    power = pump["power_fixed"] * z + pump["power_per_unit_flow"] * q * z
    return power * pump["energy_price"] * WM.ref(wm, :time_step)
end


function calc_simulation_cost(wm::AbstractCQModel)::Float64
    return max(0.0, sum(calc_pump_cost.(Ref(wm), WM.ids(wm, :pump))))
end


mutable struct Arc
    type::Symbol
    index::Int
end


function collect_arcs(wm::AbstractCQModel)::Array{Arc, 1}
    arcs = Array{Arc, 1}([])

    for symbol in [:des_pipe, :pipe, :pump, :regulator, :short_pipe, :valve]
        append!(arcs, [Arc(symbol, i) for i in collect(WM.ids(wm, symbol))])
    end

    return arcs
end


function get_pipe_flow_expression(wm::AbstractCQModel, arc::Arc)::Float64
    pipe = WM.ref(wm, :pipe, arc.index)
    q = JuMP.value(WM.var(wm, :q_pipe, arc.index))
    exponent, L = WM.ref(wm, :alpha), pipe["length"]
    head_loss, viscosity = wm.data["head_loss"], wm.data["viscosity"]
    r = WM._calc_pipe_resistance(pipe, head_loss, viscosity, 1.0, 1.0)
    return L * r * q * abs(q)^(exponent - 1.0)
end


function get_pump_flow_expression(wm::AbstractCQModel, arc::Arc)::Float64
    pump = WM.ref(wm, :pump, arc.index)
    q = max(0.0, JuMP.value(WM.var(wm, :q_pump, arc.index)))
    z = max(0.0, JuMP.value(WM.var(wm, :z_pump, arc.index)))
    c = WM._calc_head_curve_coefficients(pump)
    return -z * (c[1] * q^2 + c[2] * q + c[3])
end


function get_arc_flow_expression(wm::AbstractCQModel, arc::Arc)::Float64
    if arc.type == :pipe
        return get_pipe_flow_expression(wm, arc)
    elseif arc.type == :pump
        return get_pump_flow_expression(wm, arc)
    else
        return 0.0
    end
end


function build_right_hand_side(wm::AbstractCQModel)::Array{Float64, 1}
    node_ids, arcs = sort(collect(WM.ids(wm, :node))), collect_arcs(wm)
    reservoir_ids = sort(collect(WM.ids(wm, :reservoir)))
    tank_ids = sort(collect(WM.ids(wm, :tank)))
    node_map = Dict{Int, Int}(node_ids[i] => i for i in 1:length(node_ids))
    vector = zeros(Float64, length(arcs) + length(reservoir_ids) + length(tank_ids))

    for (i, arc) in enumerate(arcs)
        vector[i] = get_arc_flow_expression(wm, arc)
    end

    for (i, reservoir_id) in enumerate(reservoir_ids)
        node = WM.ref(wm, :node, WM.ref(wm, :reservoir, reservoir_id)["node"])
        vector[length(arcs) + i] = node["head_nominal"]
    end

    for (i, tank_id) in enumerate(tank_ids)
        tank = WM.ref(wm, :tank, tank_id)
        node = WM.ref(wm, :node, tank["node"])
        head = node["elevation"] + tank["init_level"]
        vector[length(arcs) + length(reservoir_ids) + i] = head
    end

    return vector
end


function build_incidence_matrix(wm::AbstractCQModel)::Array{Float64, 2}
    node_ids, arcs = sort(collect(WM.ids(wm, :node))), collect_arcs(wm)
    reservoir_ids = sort(collect(WM.ids(wm, :reservoir)))
    tank_ids = sort(collect(WM.ids(wm, :tank)))

    node_map = Dict{Int, Int}(node_ids[i] => i for i in 1:length(node_ids))
    num_rows = length(arcs) + length(reservoir_ids) + length(tank_ids)
    matrix = zeros(Float64, num_rows, length(node_ids))

    for (i, arc) in enumerate(arcs)
        col_fr = node_map[WM.ref(wm, arc.type, arc.index, "node_fr")]
        col_to = node_map[WM.ref(wm, arc.type, arc.index, "node_to")]

        if arc.type in [:des_pipe, :pump, :valve, :regulator]
            z = JuMP.value(wm.model[Symbol("z_" * string(arc.type))][arc.index])
            matrix[i, col_fr], matrix[i, col_to] = z, -z
        else
            matrix[i, col_fr], matrix[i, col_to] = 1.0, -1.0
        end
    end

    for (i, reservoir_id) in enumerate(reservoir_ids)
        node = WM.ref(wm, :reservoir, reservoir_id, "node")
        matrix[length(arcs) + i, node_map[node]] = 1.0
    end

    for (i, tank_id) in enumerate(tank_ids)
        node = WM.ref(wm, :tank, tank_id, "node")
        matrix[length(arcs) + length(reservoir_ids) + i, node_map[node]] = 1.0
    end

    return matrix
end


function compute_heads(incidence::Array{Float64, 2}, rhs::Array{Float64, 1})
    return incidence \ rhs
end


function flows_are_feasible(wm::AbstractCQModel)::Bool
    for symbol in [:des_pipe, :pipe, :pump, :regulator, :short_pipe, :valve]
        component_ids = sort(collect(WM.ids(wm, symbol)))
        qp = max.(0.0, JuMP.value.(WM.var.(Ref(wm), Symbol("qp_" * string(symbol)), component_ids)))
        qn = max.(0.0, JuMP.value.(WM.var.(Ref(wm), Symbol("qn_" * string(symbol)), component_ids)))

        if symbol in [:des_pipe, :pump, :regulator, :valve]
            z = JuMP.value.(WM.var.(Ref(wm), Symbol("z_" * string(symbol)), component_ids))
            vals = z .* (qp .- qn)
            lbs = [WM.ref(wm, symbol, i, "flow_min") - 1.0e-6 for i in component_ids]
            ubs = [WM.ref(wm, symbol, i, "flow_max") + 1.0e-6 for i in component_ids]
        else
            vals = qp .- qn
            lbs = [WM.ref(wm, symbol, i, "flow_min") - 1.0e-6 for i in component_ids]
            ubs = [WM.ref(wm, symbol, i, "flow_max") + 1.0e-6 for i in component_ids]
        end

        lbs_satisfied, ubs_satisfied = all(vals .>= lbs), all(vals .<= ubs)
        !(lbs_satisfied && ubs_satisfied) && return false
    end

    return true
end


function heads_are_feasible(wm::AbstractCQModel, heads::Array{Float64, 1})
    for (row_index, node_id) in enumerate(sort(collect(WM.ids(wm, :node))))
        head_min = WM.ref(wm, :node, node_id, "head_min")
        head_max = WM.ref(wm, :node, node_id, "head_max")
        head_min_satisfied = heads[row_index] >= head_min - 1.0e-6
        head_max_satisfied = heads[row_index] <= head_max + 1.0e-6
        !(head_min_satisfied && head_max_satisfied) && return false
    end

    return true
end


function check_physical_feasibility(wm::AbstractCQModel)::Bool
    if flows_are_feasible(wm)
        matrix = build_incidence_matrix(wm)
        vector = build_right_hand_side(wm)
        heads = compute_heads(matrix, vector)
        return heads_are_feasible(wm, heads)
    else
        return false
    end
end


function check_simulation_feasibility(wm::AbstractCQModel)::Bool
    solved = JuMP.termination_status(wm.model) === WM._MOI.LOCALLY_SOLVED
    almost_solved = JuMP.termination_status(wm.model) === WM._MOI.ALMOST_LOCALLY_SOLVED
    return solved || almost_solved ? check_physical_feasibility(wm) : false
end


function build_simulation_result(wm::AbstractCQModel)::SimulationResult
    q_tank = calc_simulation_tank_flows(wm)
    cost = calc_simulation_cost(wm)
    feasible = check_simulation_feasibility(wm)
    return SimulationResult(feasible, q_tank, cost)
end


function set_tank_heads!(wm::AbstractCQModel)
    for (i, tank) in WM.ref(wm, :tank)
        node = WM.ref(wm, :node, tank["node"])
        head = node["elevation"] + tank["init_level"]
        JuMP.set_value(wm.model[:h_tank][node["index"]], head)
    end
end


function set_reservoir_heads!(wm::AbstractCQModel)
    for (i, reservoir) in WM.ref(wm, :reservoir)
        node = WM.ref(wm, :node, reservoir["node"])
        JuMP.set_value(wm.model[:h_reservoir][node["index"]], node["head_nominal"])
    end
end


function set_demands!(wm::AbstractCQModel)
    for (i, demand) in WM.ref(wm, :demand)
        node = WM.ref(wm, :node, demand["node"])
        var = wm.model[:fixed_demands][node["index"]]
        JuMP.set_value(var, demand["flow_nominal"])
    end
end


function postprocess_simulation_results!(wm::AbstractCQModel)
    for (i, pump) in WM.ref(wm, :pump)
        qp = JuMP.value(WM.var(wm, :qp_pump, i))
        qn = JuMP.value(WM.var(wm, :qn_pump, i))
        z = JuMP.value(wm.model[:z_pump][i])
        JuMP.set_start_value(WM.var(wm, :qp_pump, i), qp * z)
        JuMP.set_start_value(WM.var(wm, :qn_pump, i), qn * z)
    end

    for (i, regulator) in WM.ref(wm, :regulator)
        qp = JuMP.value(WM.var(wm, :qp_regulator, i))
        qn = JuMP.value(WM.var(wm, :qn_regulator, i))
        z = JuMP.value(wm.model[:z_regulator][i])
        JuMP.set_start_value(WM.var(wm, :qp_regulator, i), qp * z)
        JuMP.set_start_value(WM.var(wm, :qn_regulator, i), qn * z)
    end

    for (i, valve) in WM.ref(wm, :valve)
        qp = JuMP.value(WM.var(wm, :qp_valve, i))
        qn = JuMP.value(WM.var(wm, :qn_valve, i))
        z = JuMP.value(wm.model[:z_valve][i])
        JuMP.set_start_value(WM.var(wm, :qp_valve, i), qp * z)
        JuMP.set_start_value(WM.var(wm, :qn_valve, i), qn * z)
    end
end


function simulate_control_setting(wm::AbstractCQModel, control_setting::ControlSetting)::SimulationResult
    WM._IM.load_timepoint!(wm.data, control_setting.network_id)
    set_tank_heads!(wm); set_reservoir_heads!(wm); set_demands!(wm);
    set_parameters_from_control_setting!(wm, control_setting)
    JuMP.optimize!(wm.model) # Optimize the internal model.
    postprocess_simulation_results!(wm)
    return build_simulation_result(wm)
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


function _instantiate_cq_model(data::Dict{String, <:Any}, optimizer)
    wm = WM.instantiate_model(data, CQWaterModel, WM.build_wf)
    JuMP.set_optimizer(wm.model, optimizer)
    return wm # Return the WaterModels object.
end


function calc_control_settings(wm::WM.AbstractWaterModel)
end


# function calc_possible_schedules(network::Dict{String, <:Any}, mip_optimizer::Any, nlp_optimizer::Any)
#     control_settings = create_all_control_settings(network)

#     network_median = deepcopy(network)
#     _set_median_tank_time_series!(network_median)

#     # Create the multinetwork version of the network.
#     network_mn = WM.make_multinetwork(network_median)

#     # Solve a relaxation of the OWF to get an estimate of tank volume time series.
#     ext = Dict(:pipe_breakpoints => 25, :pump_breakpoints => 25)
#     wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf; ext = ext)
#     result_mip = WM.optimize_model!(wm; optimizer = nlp_optimizer, relax_integrality = true)

#     # Get all component schedules across all possible time steps.
#     schedules = create_all_schedules(wm)
#     return schedules

#     # Create copies of the network to compute everything in parallel.
#     network_tmp = [deepcopy(network_median) for i in 1:Threads.nthreads()]
#     wms = [_instantiate_cq_model(network_tmp[i]) for i in 1:Threads.nthreads()]
#     return wms

#     # time = @elapsed Threads.@threads for n in sort(collect(keys(schedules)))
#     #     # Simulate schedules and filter the solutions.
#     #     WM._IM.load_timepoint!(network_tmp[Threads.threadid()], n)
#     #     results = simulate_schedules(network_tmp[Threads.threadid()], schedules[n], nlp_optimizer)
#     #     ids = filter(k -> results[k][1] === FEASIBLE_POINT, 1:length(results))
#     #     schedules[n] = (schedules[n][1], [(schedules[n][2][i], results[i][2]) for i in ids])
#     # end

#     # # Return schedules with relevant solution properties.
#     # return schedules
# end


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
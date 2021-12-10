import Random

function get_controllable_variable_indices(wm::WM.AbstractWaterModel, n::Int)
    var_ids = Array{WM._VariableIndex, 1}([])

    for comp_type in [:des_pipe, :pump, :valve, :regulator]
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
        var = WM.var(wm, symbol)[id]
        val = round(control_setting.vals[i])
        JuMP.set_value(var, val)
    end
end


function calc_simulation_tank_flows(wm::AbstractCQModel)::Dict{Int,Float64}
    return Dict{Int, Float64}(
        i => JuMP.value(WM.var(wm, :q_tank, i)) for i in WM.ids(wm, :tank))
end


function calc_pump_cost(wm::AbstractCQModel, pump_index::Int64)::Float64
    pump = WM.ref(wm, :pump, pump_index)
    q = max(0.0, JuMP.value(WM.var(wm, :q_pump, pump_index)))
    z = round(JuMP.value(wm.model[:z_pump][pump_index]))
    power = pump["power_fixed"] * z + pump["power_per_unit_flow"] * q * z
    return power * pump["energy_price"] * WM.ref(wm, :time_step)
end


function calc_simulation_cost(wm::AbstractCQModel)::Float64
    return sum(calc_pump_cost.(Ref(wm), WM.ids(wm, :pump)))
end


mutable struct Arc
    type::Symbol
    index::Int
end


function collect_arcs(wm::AbstractCQModel)::Vector{Arc}
    arcs = Vector{Arc}([])

    for symbol in [:des_pipe, :pipe, :pump, :regulator, :short_pipe, :valve]
        append!(arcs, [Arc(symbol, i) for i in sort(collect(WM.ids(wm, symbol)))])
    end

    return arcs
end


function get_pipe_flow_expression(wm::AbstractCQModel, arc::Arc)::Float64
    pipe = WM.ref(wm, :pipe, arc.index)
    q = JuMP.value(WM.var(wm, :q_pipe, arc.index))

    exponent = WM._get_exponent_from_head_loss_form(
        wm.ref[:it][WM.wm_it_sym][:head_loss])
    L = WM.ref(wm, :pipe, arc.index, "length")

    data = WM.get_wm_data(wm.data)

    head_loss = data["head_loss"]
    viscosity = data["viscosity"]
    
    base_length = get(data, "base_length", 1.0)
    base_mass = get(data, "base_mass", 1.0)
    base_time = get(data, "base_time", 1.0)

    r = WM._calc_pipe_resistance(pipe, head_loss,
        viscosity, base_length, base_mass, base_time)
        
    return L * r * sign(q) * abs(q)^exponent
end


function get_pump_flow_expression(wm::AbstractCQModel, arc::Arc)::Float64
    q = max(0.0, JuMP.value(WM.var(wm, :qp_pump, arc.index)))
    z = round(JuMP.value(WM.var(wm, :z_pump, arc.index)))
    head_curve_func = WM.ref(wm, :pump, arc.index, "head_curve_function")
    return -z * head_curve_func(q)
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


function build_right_hand_side(wm::AbstractCQModel)::Vector{Float64}
    node_ids, arcs = sort(collect(WM.ids(wm, :node))), collect_arcs(wm)
    reservoir_ids = sort(collect(WM.ids(wm, :reservoir)))
    tank_ids = sort(collect(WM.ids(wm, :tank)))
    node_map = Dict{Int, Int}(node_ids[i] => i for i in 1:length(node_ids))
    vector = zeros(Float64, length(arcs) + length(reservoir_ids) + length(tank_ids))

    for (i, arc) in enumerate(arcs)
        vector[i] = get_arc_flow_expression(wm, arc)
    end

    for (i, reservoir_id) in enumerate(reservoir_ids)
        node_id = WM.ref(wm, :reservoir, reservoir_id, "node")
        node = WM.ref(wm, :node, node_id)
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


function check_solution(wm::AbstractCQModel, rhs::Vector{Float64}, heads::Vector{Float64})::Bool
    node_ids = sort(collect(WM.ids(wm, :node)))
    node_map = Dict{Int, Int}(node_ids[i] => i for i in 1:length(node_ids))

    for (arc_index, arc) in enumerate(collect_arcs(wm))
        col_fr = node_map[WM.ref(wm, arc.type, arc.index, "node_fr")]
        col_to = node_map[WM.ref(wm, arc.type, arc.index, "node_to")]

        if arc.type in [:des_pipe, :pump, :valve, :regulator]
            z = JuMP.value(wm.model[Symbol("z_" * string(arc.type))][arc.index])
            head_difference = round(z) * (heads[col_fr] - heads[col_to])
            @assert isapprox(rhs[arc_index], head_difference; atol = 1.0e-4)
        else
            head_difference = heads[col_fr] - heads[col_to]
            @assert isapprox(rhs[arc_index], head_difference; atol = 1.0e-4)
        end
    end
    
    return true
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
            matrix[i, col_fr], matrix[i, col_to] = round(z), -round(z)
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


function compute_heads(incidence::Array{Float64, 2}, rhs::Vector{Float64})
    return LinearAlgebra.pinv(incidence) * rhs
end


function flows_are_feasible(wm::AbstractCQModel)::Bool
    for symbol in [:des_pipe, :pipe, :pump, :regulator, :short_pipe, :valve]
        component_ids = sort(collect(WM.ids(wm, symbol)))
        
        qp = max.(0.0, JuMP.value.(WM.var.(Ref(wm), Symbol("qp_" * string(symbol)), component_ids)))
        qn = max.(0.0, JuMP.value.(WM.var.(Ref(wm), Symbol("qn_" * string(symbol)), component_ids)))

        if symbol in [:des_pipe, :regulator, :valve]
            z = JuMP.value.(WM.var.(Ref(wm), Symbol("z_" * string(symbol)), component_ids))
            vals = round.(z) .* (qp .- qn)
            lbs = [WM.ref(wm, symbol, i, "flow_min") for i in component_ids]
            ubs = [WM.ref(wm, symbol, i, "flow_max") for i in component_ids]
        elseif symbol == :pump
            z = JuMP.value.(WM.var.(Ref(wm), Symbol("z_" * string(symbol)), component_ids))
            vals = round.(z) .* qp
            lbs = round.(z) .* [WM.ref(wm, symbol, i, "flow_min_forward") for i in component_ids]
            ubs = round.(z) .* [WM.ref(wm, symbol, i, "flow_max") for i in component_ids]
        else
            vals = qp .- qn
            lbs = [WM.ref(wm, symbol, i, "flow_min") for i in component_ids]
            ubs = [WM.ref(wm, symbol, i, "flow_max") for i in component_ids]
        end

        lbs_satisfied = all(vals .>= lbs .- 1.0e-6)
        ubs_satisfied = all(vals .<= ubs .+ 1.0e-6)

        if !lbs_satisfied || !ubs_satisfied
            return false
        end
    end

    return true
end


function heads_are_feasible(wm::AbstractCQModel, heads::Vector{Float64})
    for (row_index, node_id) in enumerate(sort(collect(WM.ids(wm, :node))))
        head_min = WM.ref(wm, :node, node_id, "head_min")
        head_min_satisfied = heads[row_index] >= head_min - 1.0e-6

        head_max = WM.ref(wm, :node, node_id, "head_max")
        head_max_satisfied = heads[row_index] <= head_max + 1.0e-6

        if !head_min_satisfied || !head_max_satisfied
            return false
        end
    end

    return true
end


function check_physical_feasibility(wm::AbstractCQModel)::Bool
    if flows_are_feasible(wm)
        matrix = build_incidence_matrix(wm)
        vector = build_right_hand_side(wm)
        heads = compute_heads(matrix, vector)
        valid_solution = check_solution(wm, vector, heads)
        return heads_are_feasible(wm, heads)
    else
        return false
    end
end


function check_simulation_feasibility(wm::AbstractCQModel)::Bool
    solved = JuMP.termination_status(wm.model) === WM.JuMP.MOI.LOCALLY_SOLVED
    almost_solved = JuMP.termination_status(wm.model) === WM.JuMP.MOI.ALMOST_LOCALLY_SOLVED
    return solved || almost_solved ? check_physical_feasibility(wm) : false
end


function build_simulation_result(wm::AbstractCQModel)::SimulationResult
    q_tank = calc_simulation_tank_flows(wm)
    cost = calc_simulation_cost(wm)
    feasible = check_simulation_feasibility(wm)
    return SimulationResult(feasible, q_tank, cost)
end


function set_tank_heads!(wm::AbstractCQModel)
    for (_, tank) in WM.ref(wm, :tank)
        node = WM.ref(wm, :node, tank["node"])
        head = node["elevation"] + tank["init_level"]
        var = wm.model[:h_tank][node["index"]]
        JuMP.set_value(var, head)
    end
end


function set_reservoir_heads!(wm::AbstractCQModel)
    for (_, reservoir) in WM.ref(wm, :reservoir)
        node = WM.ref(wm, :node, reservoir["node"])
        var = wm.model[:h_reservoir][node["index"]]
        JuMP.set_value(var, node["head_nominal"])
    end
end


function set_demands!(wm::AbstractCQModel)
    for (_, demand) in WM.ref(wm, :demand)
        node = WM.ref(wm, :node, demand["node"])
        var = wm.model[:fixed_demands][node["index"]]
        JuMP.set_value(var, demand["flow_nominal"])
    end
end


function set_warm_start_loose!(wm::AbstractCQModel)
    vars = JuMP.all_variables(wm.model)
    JuMP.set_start_value.(vars, 1.0e-6)
end


function simulate_control_setting(wm::AbstractCQModel, control_setting::ControlSetting)::SimulationResult
    wm_data = WM.get_wm_data(wm.data)
    WM._IM.load_timepoint!(wm_data, control_setting.network_id)

    set_tank_heads!(wm)
    set_reservoir_heads!(wm)
    set_demands!(wm)
    
    set_parameters_from_control_setting!(wm, control_setting)
    JuMP.optimize!(wm.model)
    return build_simulation_result(wm)
end


function update_data_from_simulation_result(wm::AbstractCQModel, result::SimulationResult, nw::Int)
    if !haskey(wm.data, "time_series")
        wm.data["time_series"] = Dict{String, Any}()
    end

    if !haskey(wm.data["time_series"], "tank")
        wm.data["time_series"]["tank"] = Dict{String, Any}()

        for (i, tank) in wm.data["tank"]
            wm.data["time_series"]["tank"][i] = Dict{String, Any}()
            level_array = ones(wm.data["time_series"]["num_steps"]) * tank["init_level"]
            wm.data["time_series"]["tank"][i]["init_level"] = level_array
        end
    end

    for (i, tank) in wm.data["time_series"]["tank"]
        surface_area = 0.25 * pi * wm.data["tank"][i]["diameter"]^2
        tank["init_level"][nw+1] = tank["init_level"][nw] -
            result.q_tank[parse(Int, i)] / surface_area *
            wm.data["time_series"]["time_step"]
    end
end


function simulate_control_settings_sequential(wm::AbstractCQModel, control_settings::Vector{ControlSetting})
    simulation_results = Vector{SimulationResult}([])

    for nw_id in sort(collect(keys(control_settings)))
        simulation_result = simulate_control_setting(wm, control_settings[nw_id])
        update_data_from_simulation_result(wm, simulation_result, nw_id)
        push!(simulation_results, simulation_result)

        if !simulation_result.feasible
            return simulation_results
        end
    end
        
    # Check if the tank level constraints have been violated.
    for (i, tank) in wm.data["time_series"]["tank"]
        if tank["init_level"][end] < tank["init_level"][1]
            simulation_results[end].feasible = false
        elseif tank["init_level"][end] > wm.data["tank"][i]["max_level"]
            simulation_results[end].feasible = false
        end
    end

    return simulation_results
end


function _instantiate_cq_model(data::Dict{String, <:Any}, optimizer)
    wm = WM.instantiate_model(data, CQWaterModel, WM.build_wf)
    JuMP.set_optimizer(wm.model, optimizer)
    return wm # Return the WaterModels object.
end

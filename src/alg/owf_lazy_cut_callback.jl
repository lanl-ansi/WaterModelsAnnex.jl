function _update_control_time_series!(wm::WM.AbstractWaterModel, cb_data, network)
    if !haskey(data["time_series"], "pump")
        data["time_series"]["pump"] = Dict{String, Any}()
    end

    if !haskey(data["time_series"], "regulator")
        data["time_series"]["regulator"] = Dict{String, Any}()
    end

    if !haskey(data["time_series"], "valve")
        data["time_series"]["valve"] = Dict{String, Any}()
    end

    for (k, var_id) in enumerate(schedule[1][1])
        comp_type = string(var_id.component_type)
        comp_id = string(var_id.component_index)
        comp = data[comp_type][comp_id]

        if !haskey(data["time_series"][comp_type], comp_id)
            data["time_series"][comp_type][comp_id] = Dict{String, Any}()
        end

        statuses = [Int(result[n][k]) for n in sort(collect(keys(schedule)))]
        data["time_series"][comp_type][comp_id]["status"] = Array{Int64, 1}(statuses)
        data["time_series"][comp_type][comp_id]["z_min"] = Array{Int64, 1}(statuses)
        data["time_series"][comp_type][comp_id]["z_max"] = Array{Int64, 1}(statuses)
    end
end


function set_component_statuses_from_callback_data(wm::WM.AbstractWaterModel, cb_data, network, comp_type::Symbol)
    if !haskey(network["time_series"], string(comp_type))
        network["time_series"][string(comp_type)] = Dict{String, Any}()
    end

    for (i, component) in WM.ref(wm, 1, comp_type)
        if !haskey(network["time_series"][string(comp_type)], string(i))
            status_dict = Dict{String, Any}("status" => Array{WM.STATUS, 1}([]))
            network["time_series"][string(comp_type)][string(i)] = status_dict
        else
            ts = network["time_series"][string(comp_type)][string(i)]
            network["time_series"][string(comp_type)][string(i)] =
                merge(ts, Dict("status" => Array{WM.STATUS, 1}([])))
        end
    end

    var_symbol = Symbol("z_" * string(comp_type))

    for nw in sort(collect(WM.nw_ids(wm)))
        for (i, component) in WM.ref(wm, nw, comp_type)
            z_var = WM.var(wm, nw, var_symbol, i)
            z_val = JuMP.callback_value(cb_data, z_var)
            time_series = network["time_series"][string(comp_type)][string(i)]["status"]
            push!(time_series, WM.STATUS(Int(round(z_val))))
         end
    end
end


function set_all_component_statuses_from_callback_data(wm::WM.AbstractWaterModel, cb_data, network::Dict{String, <:Any})
    for comp_type in [:pump, :regulator, :valve]
        set_component_statuses_from_callback_data(wm, cb_data, network, comp_type)
    end
end


function tank_recovery_satisfied(data::Dict{String, <:Any}, result::Dict{String, <:Any})
    for (i, tank) in data["tank"]
        node = data["node"][string(tank["node"])]
        head_min = node["elevation"] + tank["init_level"]
        head_value = result["solution"]["node"][string(node["index"])]["h"]
        head_value < head_min && (return false)
    end

    return true
end


function simulate_from_data(data::Dict{String, <:Any}, optimizer)
    result = Dict{String, Any}("primal_status" => WM._MOI.INFEASIBLE_POINT)
    cost, first_infeasible_nw, recovered = 0.0, nothing, true

    for n in 1:data["time_series"]["num_steps"]
        # Load and fix data at the current time step.
        WM._IM.load_timepoint!(data, n)
        WM.fix_all_indicators!(data)
        result = WM.solve_wf(data, CDWaterModel,
            optimizer; relax_integrality = true)

        if !feasible_simulation_result(result)
            first_infeasible_nw = n; cost = 0.0; break;
        elseif haskey(result["solution"], "pump")
            cost += sum(x["c"] for (i, x) in result["solution"]["pump"])
        end

        if n < data["time_series"]["num_steps"]
            _update_initial_tank_heads!(data, result, n + 1)
        end
    end

    if feasible_simulation_result(result)
        if !tank_recovery_satisfied(data, result)
            first_infeasible_nw = data["time_series"]["num_steps"]
            result["primal_status"] = WM._MOI.INFEASIBLE_POINT
            recovered = false
        end
    end
    
    # Return the result.
    return result, first_infeasible_nw, cost, recovered
end


function simulate_callback_solution(wm::WM.AbstractWaterModel, cb_data, network, optimizer)
    set_all_component_statuses_from_callback_data(wm, cb_data, network)
    return simulate_from_data(network, optimizer)
end


function _get_indicator_variables(wm::WM.AbstractWaterModel)
    vars = Array{WM.JuMP.VariableRef, 1}()

    for var_sym in [:z_pump, :z_regulator, :z_valve]
        for nw_id in WM.nw_ids(wm)
            append!(vars, vcat(WM.var(wm, nw_id, var_sym)...))
        end
    end

    return vars
end


function _get_pump_indicator_variables_to_nw(wm::WM.AbstractWaterModel, nw_last::Int)
    vars = Array{WM.JuMP.VariableRef, 1}()
    network_ids = sort(collect(WM.nw_ids(wm)))[1:nw_last]

    for nw_id in network_ids
        append!(vars, vcat(WM.var(wm, nw_id, :z_pump)...))
    end

    return vars
end


function _get_indicator_variables_to_nw(wm::WM.AbstractWaterModel, nw_last::Int)
    vars = Array{WM.JuMP.VariableRef, 1}()
    network_ids = sort(collect(WM.nw_ids(wm)))[1:nw_last]

    for var_sym in [:z_pump, :z_regulator, :z_valve]
        for nw_id in network_ids
            append!(vars, vcat(WM.var(wm, nw_id, var_sym)...))
        end
    end

    return vars
end


function feasible_simulation_result(result::Dict{String, <:Any})
    feasible_statuses = [WM._MOI.FEASIBLE_POINT, WM._MOI.NEARLY_FEASIBLE_POINT]
    return result["primal_status"] in feasible_statuses
end


function add_feasibility_cut!(wm::WM.AbstractWaterModel, cb_data, nw_last::Int)    
    # Collect the current integer solution into "zero" and "one" buckets.
    vars = _get_indicator_variables_to_nw(wm, nw_last)
    zero_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 0.0, vars)
    one_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 1.0, vars)

    # If the solution is not feasible (according to a comparison with WNTR), add a no-good cut.
    con = WM.JuMP.@build_constraint(sum(zero_vars) - sum(one_vars) >= 1.0 - length(one_vars))
    WM._MOI.submit(wm.model, WM._MOI.LazyConstraint(cb_data), con)
end


function add_recovery_feasibility_cut!(wm::WM.AbstractWaterModel, cb_data)    
    # Collect the current integer solution into "zero" and "one" buckets.
    vars = _get_pump_indicator_variables_to_nw(wm, length(WM.nw_ids(wm)) - 1)
    zero_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 0.0, vars)

    # If the solution is not feasible (according to a comparison with WNTR), add a no-good cut.
    con = WM.JuMP.@build_constraint(sum(zero_vars) >= 1)
    WM._MOI.submit(wm.model, WM._MOI.LazyConstraint(cb_data), con)
end


function add_objective_cut!(wm::WM.AbstractWaterModel, cb_data, objective_value::Float64)
    # Collect the current integer solution into "zero" and "one" buckets.
    vars = _get_indicator_variables(wm) # All relevant component status variables.
    zero_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 0.0, vars)
    one_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 1.0, vars)

    # Compute the objectives of the relaxation and WNTR solution.
    objective_var = WM.JuMP.variable_by_name(wm.model, "obj_aux")
    relaxed_objective = WM.JuMP.callback_value(cb_data, objective_var)

    # Add a cut to raise the objective for the solution to the true objective.
    bin_expr = objective_value * (length(one_vars) - sum(one_vars) + sum(zero_vars))
    con = WM.JuMP.@build_constraint(objective_var >= objective_value - bin_expr)
    WM._MOI.submit(wm.model, WM._MOI.LazyConstraint(cb_data), con)
end


function get_owf_lazy_cut_callback(wm::WM.AbstractWaterModel, network, optimizer)
    return function callback_function(cb_data)
        sim_time = @elapsed result_sim, nw_infeasible, cost, recovered =
            simulate_callback_solution(wm, cb_data, network, optimizer)

        if !feasible_simulation_result(result_sim) && !recovered
            add_recovery_feasibility_cut!(wm, cb_data)
        elseif !feasible_simulation_result(result_sim) && nw_infeasible !== nothing
            add_feasibility_cut!(wm, cb_data, nw_infeasible)
        elseif feasible_simulation_result(result_sim)
            add_objective_cut!(wm, cb_data, cost)
        end
    end
end


function add_owf_lazy_cut_callback!(wm::WM.AbstractWaterModel, network, optimizer)
    callback_function = get_owf_lazy_cut_callback(wm, network, optimizer)
    WM._MOI.set(wm.model, WM._MOI.LazyConstraintCallback(), callback_function)
end

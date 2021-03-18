function set_component_statuses_from_callback_data(wm::WM.AbstractWaterModel, cb_data, comp_type::Symbol)
    var_symbol = Symbol("z_" * string(comp_type))

    for (nw, nw_ref) in WM.nws(wm)
        for (i, component) in WM.ref(wm, nw, comp_type)
            z_var = WM.var(wm, nw, var_symbol, i)
            z_val = JuMP.callback_value(cb_data, z_var)
            component["status"] = WM.STATUS(Int(round(z_val)))
         end
    end
end


function set_all_component_statuses_from_callback_data(wm::WM.AbstractWaterModel, cb_data)
    for comp_type in [:pump, :regulator, :valve]
        set_component_statuses_from_callback_data(wm, cb_data, comp_type)
    end
end


function tank_recovery_satisfied(data::Dict{String, <:Any}, result::Dict{String, <:Any})
    nw_ids = sort([parse(Int, n) for n in keys(data["nw"])])

    for (i, tank) in data["nw"][string(nw_ids[1])]["tank"]
        node = data["nw"][string(nw_ids[1])]["node"][string(tank["node"])]
        head_min = node["elevation"] + tank["init_level"]
        head_value = result["solution"]["node"][string(node["index"])]["h"]
        head_value < head_min && (return false)
    end

    return true
end


function simulate_from_data(data::Dict{String, <:Any}, optimizer)
    data["multinetwork"] = false
    result = Dict{String, Any}()
    cost, first_infeasible_nw = 0.0, nothing
    nw_ids = sort([parse(Int, n) for n in keys(data["nw"])])

    for nw_id in nw_ids
        WM._IM.update_data!(data, data["nw"][string(nw_id)])
        WM.fix_all_indicators!(data)
        result = WM.solve_wf(data, CDWaterModel, optimizer; relax_integrality = true)
        WM._IM.update_data!(data["nw"][string(nw_id)], result["solution"])
        _update_initial_tank_heads!(data, result, data["nw"][string(nw_id)]["time_step"])
        
        if !feasible_simulation_result(result)
            first_infeasible_nw = nw_id; cost = 0.0; break;
        elseif haskey(result["solution"], "pump")
            cost += sum(x["c"] for (i, x) in result["solution"]["pump"])
        end
    end

    if feasible_simulation_result(result)
        if !tank_recovery_satisfied(data, result)
            first_infeasible_nw = nw_ids[end]
            result["primal_status"] = WM._MOI.INFEASIBLE_POINT
        end
    end

    # Reset data to the first time step.
    WM._IM.update_data!(data, data["nw"]["1"])
    data["multinetwork"] = true

    # Return the result.
    return result, first_infeasible_nw, cost
end


function simulate_callback_solution(wm::WM.AbstractWaterModel, cb_data, optimizer)
    set_all_component_statuses_from_callback_data(wm, cb_data)
    return simulate_from_data(wm.data, optimizer)
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


function get_owf_lazy_cut_callback(wm::WM.AbstractWaterModel, optimizer)
    return function callback_function(cb_data)
        result_sim, nw_infeasible, cost = simulate_callback_solution(wm, cb_data, optimizer)
        
        if !feasible_simulation_result(result_sim) && nw_infeasible !== nothing
            add_feasibility_cut!(wm, cb_data, nw_infeasible)
        elseif feasible_simulation_result(result_sim)
            add_objective_cut!(wm, cb_data, 0.98 * cost)
        end
    end
end


function add_owf_lazy_cut_callback!(wm::WM.AbstractWaterModel, optimizer)
    callback_function = get_owf_lazy_cut_callback(wm, optimizer)
    WM._MOI.set(wm.model, WM._MOI.LazyConstraintCallback(), callback_function)
end
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

    for nw in sort(collect(WM.nw_ids(wm)))[1:end-1]
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

        head_tm1 = result["solution"]["node"][string(node["index"])]["h"]
        q_tm1 = result["solution"]["tank"][i]["q"]
        surface_area = 0.25 * pi * tank["diameter"]^2
        head_value = head_tm1 - q_tm1 * data["time_step"] / surface_area

        head_value < head_min && (return false)
    end

    return true
end


function simulate_from_data(data::Dict{String, <:Any}, optimizer)
    result = Dict{String, Any}("primal_status" => WM._MOI.INFEASIBLE_POINT)
    cost, first_infeasible_nw, recovered = 0.0, nothing, true

    for n in 1:data["time_series"]["num_steps"] - 1
        # Load and fix data at the current time step.
        WM._IM.load_timepoint!(data, n)
        WM.fix_all_indicators!(data)
        result = WM.solve_wf(data, CDWaterModel,
            optimizer; relax_integrality = true)

        if !feasible_simulation_result(result)
            result["primal_status"] = WM._MOI.INFEASIBLE_POINT
            result["objective"] = Inf
            first_infeasible_nw = n; cost = 0.0; break;
        elseif haskey(result["solution"], "pump")
            cost += sum(x["c"] for (i, x) in result["solution"]["pump"])
        end

        result["primal_status"] = WM._MOI.FEASIBLE_POINT
        result["objective"] = 0.0

        _update_initial_tank_heads!(data, result, n + 1)
    end

    if feasible_simulation_result(result)
        WM._IM.load_timepoint!(data, 1)

        if !tank_recovery_satisfied(data, result)
            first_infeasible_nw = data["time_series"]["num_steps"] - 1
            result["primal_status"] = WM._MOI.INFEASIBLE_POINT
            result["objective"] = Inf
            recovered = false
        else
            result["primal_status"] = WM._MOI.FEASIBLE_POINT
            result["objective"] = 0.0
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
        for nw_id in sort(collect(WM.nw_ids(wm)))[1:end-1]
            append!(vars, vcat(WM.var(wm, nw_id, var_sym)...))
        end
    end

    return vars
end


function _get_direction_variables(wm::WM.AbstractWaterModel)
    vars = Array{WM.JuMP.VariableRef, 1}()

    for var_sym in [:y_des_pipe, :y_pipe, :y_short_pipe, :y_pump, :y_regulator, :y_valve]
        for nw_id in sort(collect(WM.nw_ids(wm)))[1:end-1]
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
    #feasible_statuses = [WM._MOI.FEASIBLE_POINT, WM._MOI.NEARLY_FEASIBLE_POINT]
    status_is_feasible = result["primal_status"] === WM._MOI.FEASIBLE_POINT
    objective_is_small = result["objective"] <= 1.0e-6
    return status_is_feasible && objective_is_small
end


function add_feasibility_cut!(wm::WM.AbstractWaterModel, cb_data, nw_last::Int)    
    # Collect the current integer solution into "zero" and "one" buckets.
    vars = _get_indicator_variables_to_nw(wm, nw_last)
    zero_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) < 0.5, vars)
    one_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) >= 0.5, vars)

    # If the solution is not feasible (according to a comparison with WNTR), add a no-good cut.
    con = WM.JuMP.@build_constraint(sum(zero_vars) - sum(one_vars) >= 1.0 - length(one_vars))
    WM._MOI.submit(wm.model, WM._MOI.LazyConstraint(cb_data), con)
end


function add_recovery_feasibility_cut!(wm::WM.AbstractWaterModel, cb_data)    
    # Collect the current integer solution into "zero" and "one" buckets.
    vars = _get_pump_indicator_variables_to_nw(wm, sort(collect(WM.nw_ids(wm)))[end-1])
    zero_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 0.0, vars)

    # If the solution is not feasible (according to a comparison with simulation), add a no-good cut.
    con = WM.JuMP.@build_constraint(sum(zero_vars) >= 1)
    WM._MOI.submit(wm.model, WM._MOI.LazyConstraint(cb_data), con)
end


function add_objective_cut!(wm::WM.AbstractWaterModel, cb_data, objective_value::Float64)
    # Collect the current integer solution into "zero" and "one" buckets.
    vars = _get_indicator_variables(wm) # All relevant component status variables.
    zero_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 0.0, vars)
    one_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 1.0, vars)

    # Compute the objectives of the relaxation and WNTR solution.
    obj_func = WM.JuMP.objective_function(wm.model)

    # Add a cut to raise the objective for the solution to the true objective.
    bin_expr = objective_value * (sum(one_vars) - length(one_vars) + sum(zero_vars))
    con = WM.JuMP.@build_constraint(obj_func >= objective_value - bin_expr)
    WM._MOI.submit(wm.model, WM._MOI.LazyConstraint(cb_data), con)
end


function update_master_setting_at_nw!(wm_master, setting, cb_data, nw)
    setting.network_id = nw

    for (i, var_id) in enumerate(setting.variable_indices)
        var_id.network_index = nw
        var = WM.var(wm_master, nw, var_id.variable_symbol, var_id.component_index)
        setting.vals[i] = WM.JuMP.callback_value(cb_data, var)
    end
end


function simulate_master_solution(wm_master, wm_sim, setting, cb_data)
    total_cost = 0.0
    network_ids = sort(collect(WM.nw_ids(wm_master)))

    for nw in network_ids[1:end-1]
        update_master_setting_at_nw!(wm_master, setting, cb_data, nw)
        result = simulate_control_setting(wm_sim, setting)
        !result.feasible && return 0.0, nw
        update_tank_time_series(wm_sim.data, result, nw)
        total_cost += result.cost
    end

    if !tank_levels_recovered(wm_sim.data)
        return total_cost, network_ids[end]
    else
        return total_cost, nothing
    end
end


function get_control_settings_at_nw_cb(wm::WM.AbstractWaterModel, cb_data, nw::Int)
    pump_vids = WM._VariableIndex.(nw, :pump, :z_pump, WM.ids(wm, nw, :pump))
    valve_vids = WM._VariableIndex.(nw, :valve, :z_valve, WM.ids(wm, nw, :valve))
    regulator_vids = WM._VariableIndex.(nw, :regulator, :z_regulator, WM.ids(wm, nw, :regulator))

    vids = vcat(pump_vids, valve_vids, regulator_vids)
    vars = WM._get_variable_from_index.(Ref(wm), vids)
    vals = abs.(round.(WM.JuMP.callback_value.(Ref(cb_data), vars)))
    return ControlSetting(nw, vids, vals)
end


function fix_nc_from_cb(wm_nc::WM.AbstractWaterModel, wm::WM.AbstractWaterModel, cb_data)
    nw_ids = sort(collect(WM.nw_ids(wm)))

    for nw in nw_ids
        for (i, node) in WM.ref(wm, nw, :node)
            h_cb_var = WM.var(wm, nw, :h, i)
            h_nc_var = WM.var(wm_nc, nw, :h, i)
            h_cb_val = WM.JuMP.callback_value(cb_data, h_cb_var)
            JuMP.set_start_value(h_nc_var, h_cb_val)
        end
    end
    
    for nw in nw_ids[1:end-1]
        for (a, pump) in WM.ref(wm, nw, :pump)
            z_cb_var = WM.var(wm, nw, :z_pump, a)
            z_cb_val = round(WM.JuMP.callback_value(cb_data, z_cb_var))
            z_nc_var = WM.var(wm_nc, nw, :z_pump, a)
            JuMP.set_start_value(z_nc_var, z_cb_val)
            JuMP.fix(z_nc_var, z_cb_val; force = true)

            g_cb_var = WM.var(wm, nw, :g_pump, a)
            g_nc_var = WM.var(wm_nc, nw, :g_pump, a)
            g_cb_val = WM.JuMP.callback_value(cb_data, g_cb_var)
            JuMP.set_start_value(g_nc_var, z_cb_val * g_cb_val)

            P_cb_var = WM.var(wm, nw, :P_pump, a)
            P_nc_var = WM.var(wm_nc, nw, :P_pump, a)
            P_cb_val = WM.JuMP.callback_value(cb_data, P_cb_var)
            JuMP.set_start_value(P_nc_var, z_cb_val * P_cb_val)

            qp_cb_var = WM.var(wm, nw, :qp_pump, a)
            qp_cb_val = max(0.0, WM.JuMP.callback_value(cb_data, qp_cb_var))
            qn_cb_var = WM.var(wm, nw, :qn_pump, a)
            qn_cb_val = max(0.0, WM.JuMP.callback_value(cb_data, qn_cb_var))
            q_nc_var = WM.var(wm_nc, nw, :q_pump, a)
            JuMP.set_start_value(q_nc_var, z_cb_val * (qp_cb_val - qn_cb_val))

            if nw != nw_ids[1]
                z_switch_on_cb_var = WM.var(wm, nw, :z_switch_on_pump, a)
                z_switch_on_cb_val = round(WM.JuMP.callback_value(cb_data, z_switch_on_cb_var))
                z_switch_on_nc_var = WM.var(wm_nc, nw, :z_switch_on_pump, a)
                JuMP.set_start_value(z_switch_on_nc_var, z_switch_on_cb_val)
                JuMP.fix(z_switch_on_nc_var, z_switch_on_cb_val; force = true)

                z_switch_off_cb_var = WM.var(wm, nw, :z_switch_off_pump, a)
                z_switch_off_cb_val = round(WM.JuMP.callback_value(cb_data, z_switch_off_cb_var))
                z_switch_off_nc_var = WM.var(wm_nc, nw, :z_switch_off_pump, a)
                JuMP.set_start_value(z_switch_off_nc_var, z_switch_off_cb_val)
                JuMP.fix(z_switch_off_nc_var, z_switch_off_cb_val; force = true)
            end
        end

        for (a, valve) in WM.ref(wm, nw, :valve)
            qp_cb_var = WM.var(wm, nw, :qp_valve, a)
            qp_cb_val = max(0.0, WM.JuMP.callback_value(cb_data, qp_cb_var))
            qn_cb_var = WM.var(wm, nw, :qn_valve, a)
            qn_cb_val = max(0.0, WM.JuMP.callback_value(cb_data, qn_cb_var))
            q_nc_var = WM.var(wm_nc, nw, :q_valve, a)
            JuMP.set_start_value(q_nc_var, qp_cb_val - qn_cb_val)

            z_cb_var = WM.var(wm, nw, :z_valve, a)
            z_cb_val = round(WM.JuMP.callback_value(cb_data, z_cb_var))
            z_nc_var = WM.var(wm_nc, nw, :z_valve, a)
            JuMP.set_start_value(z_nc_var, z_cb_val)
            JuMP.fix(z_nc_var, z_cb_val; force = true)
        end

        for (a, regulator) in WM.ref(wm, nw, :regulator)
            qp_cb_var = WM.var(wm, nw, :qp_regulator, a)
            qp_cb_val = max(0.0, WM.JuMP.callback_value(cb_data, qp_cb_var))
            qn_cb_var = WM.var(wm, nw, :qn_regulator, a)
            qn_cb_val = max(0.0, WM.JuMP.callback_value(cb_data, qn_cb_var))
            q_nc_var = WM.var(wm_nc, nw, :q_regulator, a)
            JuMP.set_start_value(q_nc_var, qp_cb_val - qn_cb_val)

            z_cb_var = WM.var(wm, nw, :z_regulator, a)
            z_cb_val = round(WM.JuMP.callback_value(cb_data, z_cb_var))
            z_nc_var = WM.var(wm_nc, nw, :z_regulator, a)
            JuMP.set_start_value(z_nc_var, z_cb_val)
            JuMP.fix(z_nc_var, z_cb_val; force = true)
        end

        for (a, pipe) in WM.ref(wm, nw, :pipe)
            qp_cb_var = WM.var(wm, nw, :qp_pipe, a)
            qp_cb_val = max(0.0, WM.JuMP.callback_value(cb_data, qp_cb_var))
            qn_cb_var = WM.var(wm, nw, :qn_pipe, a)
            qn_cb_val = max(0.0, WM.JuMP.callback_value(cb_data, qn_cb_var))
            q_nc_var = WM.var(wm_nc, nw, :q_pipe, a)
            start_val = abs(qp_cb_val - qn_cb_val) > 1.0e-6 ? qp_cb_val - qn_cb_val : 1.0e-6
            JuMP.set_start_value(q_nc_var, start_val)
        end

        for (a, short_pipe) in WM.ref(wm, nw, :short_pipe)
            qp_cb_var = WM.var(wm, nw, :qp_short_pipe, a)
            qp_cb_val = max(0.0, WM.JuMP.callback_value(cb_data, qp_cb_var))
            qn_cb_var = WM.var(wm, nw, :qn_short_pipe, a)
            qn_cb_val = max(0.0, WM.JuMP.callback_value(cb_data, qn_cb_var))
            q_nc_var = WM.var(wm_nc, nw, :q_short_pipe, a)
            JuMP.set_start_value(q_nc_var, qp_cb_val - qn_cb_val)
        end

        for (a, tank) in WM.ref(wm, nw, :tank)
            q_cb_var = WM.var(wm, nw, :q_tank, a)
            q_cb_val = WM.JuMP.callback_value(cb_data, q_cb_var)
            q_nc_var = WM.var(wm_nc, nw, :q_tank, a)
            JuMP.set_start_value(q_nc_var, q_cb_val)
        end

        for (a, reservoir) in WM.ref(wm, nw, :reservoir)
            q_cb_var = WM.var(wm, nw, :q_reservoir, a)
            q_cb_val = WM.JuMP.callback_value(cb_data, q_cb_var)
            q_nc_var = WM.var(wm_nc, nw, :q_reservoir, a)
            JuMP.set_start_value(q_nc_var, q_cb_val)
        end
    end
end


function simulate_callback_nc(wm_nc::WM.AbstractWaterModel, wm::WM.AbstractWaterModel, cb_data)
    fix_nc_from_cb(wm_nc, wm, cb_data)
    JuMP.optimize!(wm_nc.model)
    termination_status = JuMP.termination_status(wm_nc.model)
    cost = JuMP.objective_value(wm_nc.model)
    return termination_status, cost
end


function get_owf_lazy_cut_callback(wm::WM.AbstractWaterModel, network, optimizer, stats)
    wm_sim = _instantiate_cq_model(network, optimizer)
    network_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]

    return function callback_function(cb_data)
        control_settings = get_control_settings_at_nw_cb.(Ref(wm), cb_data, network_ids)
        stats.time_elapsed += @elapsed simulation_results =
            simulate_control_settings_sequential(wm_sim, control_settings)

        if any(x -> !x.feasible, simulation_results)
            id_infeasible = findfirst(x -> !x.feasible, simulation_results)
            nw_infeasible = network_ids[id_infeasible]
            stats.time_elapsed += @elapsed add_feasibility_cut!(wm, cb_data, nw_infeasible)
            # WM.Memento.info(LOGGER, "Infeasible solution found at step $(nw_infeasible).")
        else
            cost = sum(x.cost for x in simulation_results)
            WM.Memento.info(LOGGER, "Found feasible solution with cost $(cost).")
            stats.best_cost = cost < stats.best_cost ? cost : stats.best_cost
            # stats.time_elapsed += @elapsed add_objective_cut!(wm, cb_data, cost)
        end

        stats.num_calls += 1
    end
end


mutable struct CallbackStats
    time_elapsed::Float64
    num_calls::Int64
    best_cost::Float64
end


function add_owf_lazy_cut_callback!(wm::WM.AbstractWaterModel, network, optimizer)
    callback_stats = CallbackStats(0.0, 0, Inf)
    callback_function = get_owf_lazy_cut_callback(wm, network, optimizer, callback_stats)
    WM._MOI.set(wm.model, WM._MOI.LazyConstraintCallback(), callback_function)
    return callback_stats
end

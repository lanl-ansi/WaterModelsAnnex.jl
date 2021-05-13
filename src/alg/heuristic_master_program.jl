function heuristic_master_program_indicator_variables!(model::JuMP.Model, settings)
    network_ids = sort(unique([x.network_id for x in settings]))
    return Dict{Int, Any}(nw => JuMP.@variable(model,
        [i in 1:length(settings[1].vals)],
        binary = true) for nw in network_ids)
end


function add_master_program_symmetry_constraints!(wm::WM.AbstractWaterModel, model::JuMP.Model, z, settings)
    network_ids = sort(unique([x.network_id for x in settings]))
    setting_vars = settings[1].variable_indices

    for nw in network_ids
        for (pump_group_id, pump_group) in WM.ref(wm, nw, :pump_group)
            pump_ids = sort(collect(pump_group["pump_indices"]))
            pump_var_id_ids = [findfirst(x -> x.variable_symbol == :z_pump &&
                x.component_index == i, setting_vars) for i in pump_ids]
            
             for i in 1:length(pump_var_id_ids) - 1
                var_1 = z[nw][pump_var_id_ids[i]]
                var_2 = z[nw][pump_var_id_ids[i+1]]
                JuMP.@constraint(model, var_1 >= var_2)
             end
        end
    end
end


function add_master_program_switch_constraints!(wm::WM.AbstractWaterModel, model::JuMP.Model, z, x, settings)
    network_ids = sort(unique([x.network_id for x in settings]))
    setting_vars = settings[1].variable_indices

    # for nw in network_ids
    #     for (pump_group_id, pump_group) in WM.ref(wm, nw, :pump_group)
    #         pump_ids = sort(collect(pump_group["pump_indices"]))
    #         pump_var_id_ids = [findfirst(x -> x.variable_symbol == :z_pump &&
    #             x.component_index == i, setting_vars) for i in pump_ids]
            
    #          for i in 1:length(pump_var_id_ids) - 1
    #             var_1 = z[nw][pump_var_id_ids[i]]
    #             var_2 = z[nw][pump_var_id_ids[i+1]]
    #             JuMP.@constraint(model, var_1 >= var_2)
    #          end
    #     end
    # end
end


function heuristic_master_program_delta_parameters(settings, weights)
    network_ids = unique([x.network_id for x in settings])
    delta = Dict{Int, Any}(nw => zeros(length(settings[1].vals)) for nw in network_ids)

    for (i, setting) in enumerate(settings)
        delta[setting.network_id] += weights[i] * setting.vals
    end

    return delta
end


function choose_beta_val(x::Float64; randomize::Bool = true)
    if isapprox(x, 0.0; atol = 1.0e-4)
        return 0.0
    elseif isapprox(x, 1.0; atol = 1.0e-4)
        return 1.0
    elseif randomize
        return Random.shuffle([0.0, x, 1.0])[1]
    else
        return x
    end
end


function heuristic_master_program_beta_parameters(delta; randomize::Bool = true)
    Beta = deepcopy(delta)

    for (nw, entry) in Beta
        Beta[nw] = [choose_beta_val(x; randomize = randomize) for x in entry]
    end

    return Beta
end


function heuristic_master_program_objective!(model, z, Beta)
    network_indices = sort(collect(keys(z)))
    cost_1 = sum(sum((Beta[n][i] - z[n][i])^2 for i in 1:length(z[n])) for n in network_indices)
    cost_2 = sum((sum(Beta[n][i] for n in network_indices) -
        sum(z[n][i] for n in network_indices))^2 for i in 1:length(z[1]))
    return JuMP.@objective(model, WM._MOI.MIN_SENSE, cost_1 + cost_2)
end


function optimize_heuristic_master_program!(model::JuMP.Model, z)
    JuMP.optimize!(model)
    return Dict{Int, Any}(n => JuMP.value.(z[n]) for n in keys(z))
end


function add_master_program_lazy_callback!(network, model::JuMP.Model, z, settings, nlp_optimizer)
    callback = get_master_program_lazy_callback(network, model, z, settings, nlp_optimizer)
    WM._MOI.set(model, WM._MOI.LazyConstraintCallback(), callback)
end


function update_setting_at_nw!(setting, cb_data, z, nw)
    setting.network_id = nw
    setting.vals = WM.JuMP.callback_value.(Ref(cb_data), z)
end


function update_setting_at_nw!(setting, z, nw)
    setting.network_id = nw
    setting.vals = WM.JuMP.value.(z)
end


function update_tank_time_series(network::Dict{String, <:Any}, result::SimulationResult, nw::Int)
    tank_ts = network["time_series"]["tank"]
    time_step = network["time_step"]

    for (tank_index, q_tank) in result.q_tank
        tank = network["tank"][string(tank_index)]
        coeff = time_step / (0.25 * pi * tank["diameter"]^2)
        tank["init_level"] -= coeff * q_tank
        tank_ts[string(tank_index)]["init_level"][nw+1] =
            tank_ts[string(tank_index)]["init_level"][nw] - coeff * q_tank
    end
end


function tank_levels_recovered(network::Dict{String, <:Any})
    for (tank_index, tank_ts) in network["time_series"]["tank"]
        tank_ts["init_level"][end] < tank_ts["init_level"][1] && return false
    end

    return true
end


function simulate_heuristic_master_solution(wm, model, cb_data, z, settings)
    network_ids = sort(collect(keys(z)))
    setting = deepcopy(settings[1])

    for nw in sort(collect(keys(z)))
        update_setting_at_nw!(setting, cb_data, z[nw], nw)
        result = simulate_control_setting(wm, setting)
        !result.feasible && return nw
        update_tank_time_series(wm.data, result, nw)
    end

    if !tank_levels_recovered(wm.data)
        return network_ids[end-1]
    else
        return nothing
    end
end


function simulate_heuristic_master_solution(wm, z, settings)
    network_ids = sort(collect(keys(z)))
    setting = deepcopy(settings[1])

    for nw in network_ids
        update_setting_at_nw!(setting, z[nw], nw)
        result = simulate_control_setting(wm, setting)
        !result.feasible && return nw
        update_tank_time_series(wm.data, result, nw)
    end

    if !tank_levels_recovered(wm.data)
        return network_ids[end]
    else
        return nothing
    end
end


function get_master_program_lazy_callback(network, model::JuMP.Model, z, settings, nlp_optimizer)
    wm = _instantiate_cq_model(network, nlp_optimizer)

    return function callback_function(cb_data)
        infeasible_nw = simulate_heuristic_master_solution(wm, model, cb_data, z, settings)

        if infeasible_nw !== nothing
            # Collect the current integer solution into "zero" and "one" buckets.
            vars = vcat([collect(z[nw]) for nw in 1:infeasible_nw]...)
            zero_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 0.0, vars)
            one_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 1.0, vars)
        
            # If the solution is not feasible (according to a simulation comparison), add a no-good cut.
            con = WM.JuMP.@build_constraint(sum(zero_vars) - sum(one_vars) >= 1.0 - length(one_vars))
            WM._MOI.submit(model, WM._MOI.LazyConstraint(cb_data), con)
        end
    end
end


function solve_heuristic_master_program(wm, network, settings, weights, optimizer, nlp_optimizer; max_iterations = 10)
    Random.seed!(0) # This can be commented out.
    model, num_iterations = JuMP.Model(optimizer), 0

    z = heuristic_master_program_indicator_variables!(model, settings)
    # x = heuristic_master_program_indicator_variables!(model, settings)
    delta = heuristic_master_program_delta_parameters(settings, weights)
    Beta = heuristic_master_program_beta_parameters(delta; randomize = false)
    objective = heuristic_master_program_objective!(model, z, Beta)
    add_master_program_symmetry_constraints!(wm, model, z, settings)
    # add_master_program_switch_constraints!(wm, model, z, x, settings)
    add_master_program_lazy_callback!(network, model, z, settings, nlp_optimizer)

    while num_iterations <= max_iterations
        JuMP.optimize!(model)
        JuMP.primal_status(model) == WM._MOI.FEASIBLE_POINT && break
        WM.Memento.info(LOGGER, "Reattempting heuristic solution discovery.")

        model = JuMP.Model(optimizer)
        z = heuristic_master_program_indicator_variables!(model, settings)
        # x = heuristic_master_program_indicator_variables!(model, settings)
        delta = heuristic_master_program_delta_parameters(settings, weights)
        Beta = heuristic_master_program_beta_parameters(delta; randomize = true)
        objective = heuristic_master_program_objective!(model, z, Beta)
        add_master_program_symmetry_constraints!(wm, model, z, settings)
        # add_master_program_switch_constraints!(wm, model, z, x, settings)
        add_master_program_lazy_callback!(network, model, z, settings, nlp_optimizer)

        num_iterations += 1
    end

    if JuMP.primal_status(model) == WM._MOI.FEASIBLE_POINT
        return Dict{Int, Any}(n => JuMP.value.(z[n]) for n in keys(z))
    else
        return nothing
    end
end

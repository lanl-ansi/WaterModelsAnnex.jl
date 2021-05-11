function update_tank_time_series_from_callback!(wm::WM.AbstractWaterModel, network, cb_data)
    for nw in sort(collect(WM.nw_ids(wm)))
        for (i, tank) in WM.ref(wm, nw, :tank)
            h_var = WM.var(wm, nw, :h, tank["node"])
            h_val = WM.JuMP.callback_value(cb_data, h_var)
            level = h_val - WM.ref(wm, nw, :node, tank["node"])["elevation"]
            network["time_series"]["tank"][string(i)]["init_level"][nw] = level
        end
    end
end


function submit_heuristic_solution_from_setting!(wm, cb_data, network, setting, nlp_optimizer)
    wm_sim = _instantiate_cq_model(network, nlp_optimizer)
    network_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]
    vars = Array{JuMP.VariableRef, 1}([])
    vals = Array{Float64, 1}([])
    total_cost = 0.0

    for nw in network_ids
        result = simulate_control_setting(wm_sim, setting[nw])
        total_cost += result.cost

        for (i, pump) in WM.ref(wm_sim, :pump)
            q_pump = JuMP.value(WM.var(wm_sim, :qp_pump, i))
            z_pump = JuMP.value(WM.var(wm_sim, :z_pump, i))
            
            push!(vars, WM.var(wm, nw, :z_pump, i))
            push!(vals, round(JuMP.value(WM.var(wm_sim, :z_pump, i))))
            push!(vars, WM.var(wm, nw, :y_pump, i))
            push!(vals, q_pump * z_pump > 0.0 ? 1.0 : 0.0)
        end

        for (i, pipe) in WM.ref(wm_sim, :pipe)
            qp_pipe = JuMP.value(WM.var(wm_sim, :qp_pipe, i))
            qn_pipe = JuMP.value(WM.var(wm_sim, :qn_pipe, i))
            push!(vars, WM.var(wm, nw, :y_pipe, i))
            push!(vals, qp_pipe - qn_pipe > 0.0 ? 1.0 : 0.0)
        end

        @assert result.feasible
        update_tank_time_series(wm_sim.data, result, nw)
    end

    accepted = WM._MOI.submit(wm.model, MOI.HeuristicSolution(cb_data), vars, vals)
    println(total_cost, " ", accepted)
end


function get_owf_heuristic_callback(wm::WM.AbstractWaterModel, network, nlp_optimizer, mip_optimizer)
    num_nodes = 0 # Initialize the number of nodes thus far explored.

    return function callback_function(cb_data)
        num_nodes += 1 # Update the number of nodes explored.

        if num_nodes % 100 == 0 # Run a heuristic every 500 nodes.
            settings = create_all_control_settings(wm)
            update_tank_time_series_from_callback!(wm, network, cb_data)
            results = simulate_control_settings(network, settings, nlp_optimizer)
            filter_feasible_control_settings!(settings, results)
            weights = solve_heuristic_linear_program(wm, settings, results, nlp_optimizer)
            control_sol = solve_heuristic_master_program(wm, network, settings,
                weights, mip_optimizer, nlp_optimizer; max_iterations = 10)

            if control_sol !== nothing
                setting = control_settings_from_solution(control_sol, settings)
                submit_heuristic_solution_from_setting!(wm, cb_data, network, setting, nlp_optimizer)
            end
        end
    end
end


function add_owf_heuristic_callback!(wm::WM.AbstractWaterModel, network, nlp_optimizer, mip_optimizer)
    callback_function = get_owf_heuristic_callback(wm, network, nlp_optimizer, mip_optimizer)
    WM._MOI.set(wm.model, WM._MOI.HeuristicCallback(), callback_function)
end
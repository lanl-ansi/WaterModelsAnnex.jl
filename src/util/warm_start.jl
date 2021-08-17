function set_warm_start_from_setting!(wm, network, settings, optimizer)
    wm_sim = _instantiate_cq_model(network, optimizer)
    node_ids = sort(collect(WM.ids(wm_sim, :node)))
    node_map = Dict(node_ids[i] => i for i in 1:length(node_ids))
    nw_end = sort(collect(WM.nw_ids(wm)))[end]
    total_cost = 0.0

    for nw in sort(collect(WM.nw_ids(wm)))
        for (i, node) in WM.ref(wm, nw, :node)
            head_mid = 0.5 * (node["head_min"] + node["head_max"])
            h_var = WM.var(wm, nw, :h, i) 
            JuMP.set_start_value(h_var, head_mid)
        end
    end

    for nw in sort([x.network_id for x in settings])
        result = simulate_control_setting(wm_sim, settings[nw])

        for k in 1:length(settings[nw].variable_indices)
            vid = settings[nw].variable_indices[k]
            val = settings[nw].vals[k]
            var = WM.var(wm, nw, vid.variable_symbol, vid.component_index)
            JuMP.set_start_value(var, val)
         end

        for comp_type in [:des_pipe, :pipe, :pump, :regulator, :short_pipe, :valve]
            for (i, comp) in WM.ref(wm_sim, comp_type)
                qp = WM.var(wm_sim, Symbol("qp_" * string(comp_type)), i)
                qn = WM.var(wm_sim, Symbol("qn_" * string(comp_type)), i)
                
                if comp_type in [:pump, :regulator, :valve]
                    z = wm_sim.model[Symbol("z_" * string(comp_type))][i]
                    qp_val = max(0.0, JuMP.value(qp) * round(JuMP.value(z)))
                    qn_val = max(0.0, JuMP.value(qn) * round(JuMP.value(z)))
                else
                    qp_val, qn_val = max(0.0, JuMP.value(qp)), max(0.0, JuMP.value(qn))
                end

                q_val = qp_val - qn_val
                y_val = q_val >= 0.0 ? 1.0 : 0.0
                qp_val, qn_val = q_val * y_val, -q_val * (1.0 - y_val)

                qp_var = WM.var(wm, nw, Symbol("qp_" * string(comp_type)), i)
                qn_var = WM.var(wm, nw, Symbol("qn_" * string(comp_type)), i)
                y_var = WM.var(wm, nw, Symbol("y_" * string(comp_type)), i)

                JuMP.set_start_value(qp_var, qp_val)
                JuMP.set_start_value(qn_var, qn_val)
                JuMP.set_start_value(y_var, y_val)
            end
        end

        matrix = build_incidence_matrix(wm_sim)
        vector = build_right_hand_side(wm_sim)
        heads = compute_heads(matrix, vector)

        for comp_type in [:pipe, :pump, :regulator, :short_pipe, :valve]
            for (a, comp) in WM.ref(wm_sim, comp_type)
                row_i, row_j = node_map[comp["node_fr"]], node_map[comp["node_to"]]
                dhp_val = max(0.0, heads[row_i] - heads[row_j])
                dhn_val = max(0.0, heads[row_j] - heads[row_i])

                if comp_type in [:pipe]
                    dhp = WM.var(wm, nw, Symbol("dhp_" * string(comp_type)), a)
                    dhn = WM.var(wm, nw, Symbol("dhn_" * string(comp_type)), a)

                    JuMP.set_start_value(dhp, dhp_val)
                    JuMP.set_start_value(dhn, dhn_val)
                elseif comp_type in [:pump]
                    z = wm_sim.model[Symbol("z_" * string(comp_type))][a]
                    g_val = dhn_val * round(JuMP.value(z))
                    g = WM.var(wm, nw, Symbol("g_" * string(comp_type)), a)
                    JuMP.set_start_value(g, g_val)
                end
            end
        end

        for node_id in WM.ids(wm, nw, :node)
            var = WM.var(wm, nw, :h, node_id)
            JuMP.set_start_value(var, heads[node_map[node_id]])
        end

        for (i, demand) in WM.ref(wm, nw, :dispatchable_demand)
            qd = WM.var(wm, nw, :q_demand, i)
            qd_sol = JuMP.value(WM.var(wm_sim, :q_demand, i))
            JuMP.set_start_value(qd, qd_sol)
        end

        for (i, reservoir) in WM.ref(wm, nw, :reservoir)
            qr = WM.var(wm, nw, :q_reservoir, i)
            qr_sol = JuMP.value(WM.var(wm_sim, :q_reservoir, i))
            JuMP.set_start_value(qr, qr_sol)
        end

        for (i, tank) in WM.ref(wm_sim, :tank)
            coeff = WM.ref(wm, nw, :time_step) / (0.25 * pi * tank["diameter"]^2)
            qt_sol = JuMP.value(WM.var(wm_sim, :q_tank, i))

            if nw + 1 <= nw_end
                wm_sim.data["time_series"]["tank"][string(i)]["init_level"][nw+1] =
                wm_sim.data["time_series"]["tank"][string(i)]["init_level"][nw] - coeff * qt_sol
                #println(wm_sim.data["time_series"]["tank"][string(i)]["init_level"])
            end

            qt = WM.var(wm, nw, :q_tank, i)
            JuMP.set_start_value(qt, qt_sol)

            if nw + 1 <= nw_end
                h_var = WM.var(wm, nw + 1, :h, tank["node"])
                h_val = heads[node_map[tank["node"]]] - coeff * qt_sol
                JuMP.set_start_value(h_var, h_val)
            end
        end

        for (i, pump) in WM.ref(wm, nw, :pump)
            q = JuMP.start_value(WM.var(wm, nw, :qp_pump, i))
            z = round(JuMP.start_value(WM.var(wm, nw, :z_pump, i)))
            P = pump["power_fixed"] * z + pump["power_per_unit_flow"] * q * z
            JuMP.set_start_value(WM.var(wm, nw, :P_pump, i), P)
        end

        total_cost += calc_simulation_cost(wm_sim)
    end

    # Start with the first network, representing the initial time step.
    network_ids = sort(collect(WM.nw_ids(wm)))
    n_1 = network_ids[1]

    # # Constraints on pump switches.
    # for n_2 in network_ids[2:end-1]
    #     # Constraints the define the pump switch variables.
    #     for (a, pump) in WM.ref(wm, :pump; nw = n_1)
    #         z_1_val = JuMP.start_value(WM.var(wm, n_1, :z_pump, a))
    #         z_2_val = JuMP.start_value(WM.var(wm, n_2, :z_pump, a))

    #         switch_on = max(0.0, z_2_val - z_1_val)
    #         JuMP.set_start_value(WM.var(wm, n_2, :z_switch_on_pump, a), switch_on)

    #         switch_off = max(0.0, z_1_val - z_2_val)
    #         JuMP.set_start_value(WM.var(wm, n_2, :z_switch_off_pump, a), switch_off)
    #     end

    #     # Update the first network used for integration.
    #     n_1 = n_2
    # end

    # Return the cost of the warm-started solution.
    return total_cost
end
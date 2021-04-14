function _populate_oa_dict(wm::WM.AbstractWaterModel, comp_type::Symbol)
    vals = Dict{Int, Any}(nw => Dict{Int, Array}() for nw in WM.nw_ids(wm))

    for nw in WM.nw_ids(wm)
        for (i, comp) in WM.ref(wm, nw, comp_type)
            vals[nw][i] = Array{Float64, 1}([])
        end
    end

    return vals
end

function get_owf_user_cut_callback(wm::WM.AbstractWaterModel)
    head_loss, viscosity = wm.data["head_loss"], wm.data["viscosity"]
    base_length, base_time = 1.0, 1.0
    exponent = WM._get_exponent_from_head_loss_form(head_loss)

    qp_pipe_vals = _populate_oa_dict(wm, :pipe)
    qn_pipe_vals = _populate_oa_dict(wm, :pipe)
    qp_pump_vals = _populate_oa_dict(wm, :pump)

    return function callback_function(cb_data)
        for nw in WM.nw_ids(wm)
            for (a, pipe) in WM.ref(wm, nw, :pipe)
                y = WM.var(wm, nw, :y_pipe, a)
                y_sol = WM.JuMP.callback_value(cb_data, y)
                r = WM._calc_pipe_resistance(pipe, head_loss, viscosity, base_length, base_time)

                if y_sol >= 0.5
                    qp = WM.var(wm, nw, :qp_pipe, a)
                    qp_sol = WM.JuMP.callback_value(cb_data, qp)

                    dhp = WM.var(wm, nw, :dhp_pipe, a)
                    dhp_sol = WM.JuMP.callback_value(cb_data, dhp)
                    dhp_pred = pipe["length"] * r * max(0.0, qp_sol)^exponent

                    if length(qp_pipe_vals[nw][a]) > 0
                        qp_diff = maximum(abs.(qp_pipe_vals[nw][a] .- qp_sol))
                    else
                        qp_diff = Inf
                    end

                    if abs(dhp_sol - dhp_pred) > 1.0e-4 && qp_diff > 1.0e-4
                        push!(qp_pipe_vals[nw][a], max(0.0, qp_sol))
                        lhs = WM._calc_head_loss_oa(qp, y, qp_sol, exponent)
                        con = JuMP.@build_constraint(r * lhs <= dhp / pipe["length"])
                        WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
                    end
                else
                    qn = WM.var(wm, nw, :qn_pipe, a)
                    qn_sol = WM.JuMP.callback_value(cb_data, qn)
                
                    dhn = WM.var(wm, nw, :dhn_pipe, a)
                    dhn_sol = WM.JuMP.callback_value(cb_data, dhn)
                    dhn_pred = pipe["length"] * r * max(0.0, qn_sol)^exponent

                    if length(qn_pipe_vals[nw][a]) > 0
                        qn_diff = maximum(abs.(qn_pipe_vals[nw][a] .- qn_sol))
                    else
                        qn_diff = Inf
                    end

                    if abs(dhn_sol - dhn_pred) > 1.0e-4 && qn_diff > 1.0e-4
                        push!(qn_pipe_vals[nw][a], max(0.0, qn_sol))
                        lhs = WM._calc_head_loss_oa(qn, 1.0 - y, qn_sol, exponent)
                        con = JuMP.@build_constraint(r * lhs <= dhn / pipe["length"])
                        WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
                    end
                end
            end

            for (a, pump) in WM.ref(wm, nw, :pump)
                z = WM.var(wm, nw, :z_pump, a)
                z_sol = WM.JuMP.callback_value(cb_data, z)

                if z_sol >= 0.5
                    qp = WM.var(wm, nw, :qp_pump, a)
                    qp_sol = WM.JuMP.callback_value(cb_data, qp)

                    g = WM.var(wm, nw, :g_pump, a)
                    g_sol = WM.JuMP.callback_value(cb_data, g)

                    # Calculate the head curve function and its derivative.
                    head_curve_func = WM._calc_head_curve_function(pump)
                    head_curve_deriv = WM._calc_head_curve_derivative(pump)

                    g_pred = head_curve_func(qp_sol)
                    df_pred = head_curve_deriv(qp_sol)

                    if length(qp_pump_vals[nw][a]) > 0
                        qp_diff = maximum(abs.(qp_pump_vals[nw][a] .- qp_sol))
                    else
                        qp_diff = Inf
                    end

                    if abs(g_sol - g_pred) > 1.0e-4 && qp_diff > 1.0e-4
                        push!(qp_pump_vals[nw][a], max(0.0, qp_sol))
                        con = JuMP.@build_constraint(g <= g_pred * z + df_pred * (qp - qp_sol * z))
                        WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
                    end
                end
            end
        end
    end
end


function get_owf_user_cut_callback(wm::AbstractLRDXModel)
    head_loss, viscosity = wm.data["head_loss"], wm.data["viscosity"]
    base_length, base_time = 1.0, 1.0
    exponent = WM._get_exponent_from_head_loss_form(head_loss)

    qp_pipe_vals = _populate_oa_dict(wm, :pipe)
    dhp_pipe_vals = _populate_oa_dict(wm, :pipe)
    qn_pipe_vals = _populate_oa_dict(wm, :pipe)
    dhn_pipe_vals = _populate_oa_dict(wm, :pipe)
    qp_pump_vals = _populate_oa_dict(wm, :pump)
    g_pump_vals = _populate_oa_dict(wm, :pump)

    return function callback_function(cb_data)
        for nw in WM.nw_ids(wm)
            for (a, pipe) in WM.ref(wm, nw, :pipe)
                y = WM.var(wm, nw, :y_pipe, a)
                y_sol = WM.JuMP.callback_value(cb_data, y)
                r = WM._calc_pipe_resistance(pipe, head_loss, viscosity, base_length, base_time)
                L_r, r_r = pipe["length"]^(-1.0 / exponent), r^(-1.0 / exponent)

                if y_sol >= 0.5
                    qp = WM.var(wm, nw, :qp_pipe, a)
                    qp_sol = WM.JuMP.callback_value(cb_data, qp)

                    qp_nl = WM.var(wm, nw, :qp_nl_pipe, a)
                    qp_nl_sol = WM.JuMP.callback_value(cb_data, qp_nl)
                    qp_nl_pred = pipe["length"] * r * 1.0 / (1.0 + exponent) * max(0.0, qp_sol)^(1.0 + exponent)

                    dhp = WM.var(wm, nw, :dhp_pipe, a)
                    dhp_sol = WM.JuMP.callback_value(cb_data, dhp)
                    dhp_pred = pipe["length"] * r * max(0.0, qp_sol)^exponent

                    if length(qp_pipe_vals[nw][a]) > 0
                        qp_diff = maximum(abs.(qp_pipe_vals[nw][a] .- qp_sol))
                    else
                        qp_diff = Inf
                    end

                    dhp_diff = abs(dhp_sol - dhp_pred) > 1.0e-4
                    qp_nl_diff = abs(qp_nl_sol - qp_nl_pred) > 1.0e-4

                    if (dhp_diff || qp_nl_diff) && qp_diff > 1.0e-4
                        push!(qp_pipe_vals[nw][a], max(0.0, qp_sol))

                        lhs = WM._calc_head_loss_oa(qp, y, qp_sol, exponent)
                        con = JuMP.@build_constraint(r * lhs <= dhp / pipe["length"])
                        WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)

                        lhs = _calc_pipe_flow_integrated_oa(qp, y, qp_sol, exponent)
                        con = JuMP.@build_constraint(r * lhs <= qp_nl / pipe["length"])
                        WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
                    end

                    dhp_nl = WM.var(wm, nw, :dhp_nl_pipe, a)
                    dhp_nl_sol = WM.JuMP.callback_value(cb_data, dhp_nl)
                    dhp_nl_pred = L_r * r_r * exponent / (1.0 + exponent) * max(0.0, dhp_sol)^(1.0 + 1.0 / exponent)

                    if length(dhp_pipe_vals[nw][a]) > 0
                        dhp_diff = maximum(abs.(dhp_pipe_vals[nw][a] .- dhp_sol))
                    else
                        dhp_diff = Inf
                    end

                    if abs(dhp_nl_sol - dhp_nl_pred) > 1.0e-4 && dhp_diff > 1.0e-4
                        push!(dhp_pipe_vals[nw][a], max(0.0, dhp_sol))
                        lhs = _calc_pipe_head_integrated_oa(dhp, y, dhp_sol, exponent)
                        con = JuMP.@build_constraint(r_r * lhs <= dhp_nl / L_r)
                        WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
                    end
                else
                    qn = WM.var(wm, nw, :qn_pipe, a)
                    qn_sol = WM.JuMP.callback_value(cb_data, qn)

                    qn_nl = WM.var(wm, nw, :qn_nl_pipe, a)
                    qn_nl_sol = WM.JuMP.callback_value(cb_data, qn_nl)
                    qn_nl_pred = pipe["length"] * r * 1.0 / (1.0 + exponent) * max(0.0, qn_sol)^(1.0 + exponent)

                    dhn = WM.var(wm, nw, :dhn_pipe, a)
                    dhn_sol = WM.JuMP.callback_value(cb_data, dhn)
                    dhn_pred = pipe["length"] * r * max(0.0, qn_sol)^exponent

                    if length(qn_pipe_vals[nw][a]) > 0
                        qn_diff = maximum(abs.(qn_pipe_vals[nw][a] .- qn_sol))
                    else
                        qn_diff = Inf
                    end

                    dhn_diff = abs(dhn_sol - dhn_pred) > 1.0e-4
                    qn_nl_diff = abs(qn_nl_sol - qn_nl_pred) > 1.0e-4

                    if (dhn_diff || qn_nl_diff) && qn_diff > 1.0e-4
                        push!(qn_pipe_vals[nw][a], max(0.0, qn_sol))

                        lhs = WM._calc_head_loss_oa(qn, 1.0 - y, qn_sol, exponent)
                        con = JuMP.@build_constraint(r * lhs <= dhn / pipe["length"])
                        WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)

                        lhs = _calc_pipe_flow_integrated_oa(qn, 1.0 - y, qn_sol, exponent)
                        con = JuMP.@build_constraint(r * lhs <= qn_nl / pipe["length"])
                        WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
                    end

                    dhn_nl = WM.var(wm, nw, :dhn_nl_pipe, a)
                    dhn_nl_sol = WM.JuMP.callback_value(cb_data, dhn_nl)
                    dhn_nl_pred = L_r * r_r * exponent / (1.0 + exponent) * max(0.0, dhn_sol)^(1.0 + 1.0 / exponent)

                    if length(dhn_pipe_vals[nw][a]) > 0
                        dhn_diff = maximum(abs.(dhn_pipe_vals[nw][a] .- dhn_sol))
                    else
                        dhn_diff = Inf
                    end

                    if abs(dhn_nl_sol - dhn_nl_pred) > 1.0e-4 && dhn_diff > 1.0e-4
                        push!(dhn_pipe_vals[nw][a], max(0.0, dhn_sol))
                        lhs = _calc_pipe_head_integrated_oa(dhn, 1.0 - y, dhn_sol, exponent)
                        con = JuMP.@build_constraint(r_r * lhs <= dhn_nl / L_r)
                        WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
                    end
                end
            end

            for (a, pump) in WM.ref(wm, nw, :pump)
                z = WM.var(wm, nw, :z_pump, a)
                z_sol = WM.JuMP.callback_value(cb_data, z)

                if z_sol >= 0.5
                    qp = WM.var(wm, nw, :qp_pump, a)
                    qp_sol = WM.JuMP.callback_value(cb_data, qp)

                    g = WM.var(wm, nw, :g_pump, a)
                    z = WM.var(wm, nw, :z_pump, a)
                    g_sol = WM.JuMP.callback_value(cb_data, g)

                    # Calculate the head curve function and its derivative.
                    head_curve_func = WM._calc_head_curve_function(pump)
                    head_curve_deriv = WM._calc_head_curve_derivative(pump)

                    g_pred = head_curve_func(qp_sol)
                    df_pred = head_curve_deriv(qp_sol)

                    if length(qp_pump_vals[nw][a]) > 0
                        qp_diff = maximum(abs.(qp_pump_vals[nw][a] .- qp_sol))
                    else
                        qp_diff = Inf
                    end

                    c = WM._calc_head_curve_coefficients(pump)
                    qp_nl = WM.var(wm, nw, :qp_nl_pump, a)
                    qp_nl_sol = WM.JuMP.callback_value(cb_data, qp_nl)
                    qp_nl_pred = (1.0 / 3.0) * c[1] * qp_sol^3 + 0.5 * c[2] * qp_sol^2 + c[3] * qp_sol

                    g_diff = abs(g_sol - g_pred) > 1.0e-4
                    qp_nl_diff = abs(qp_nl_sol - qp_nl_pred) > 1.0e-4

                    if (g_diff || qp_nl_diff) && qp_diff > 1.0e-4
                        push!(qp_pump_vals[nw][a], max(0.0, qp_sol))

                        con = JuMP.@build_constraint(g <= g_pred * z + df_pred * (qp - qp_sol * z))
                        WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)

                        rhs = _calc_pump_flow_integrated_oa(qp, z, qp_sol, c)
                        con = JuMP.@build_constraint(qp_nl <= rhs)
                        WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
                    end

                    g_nl = WM.var(wm, nw, :g_nl_pump, a)
                    g_nl_sol = WM.JuMP.callback_value(cb_data, g_nl)
                    g_nl_pred = ((-sqrt(-4.0 * c[1] * c[3] + 4.0 * c[1] * g_sol + c[2]^2) -
                        c[2])^3 / (24.0 * c[1]^2) + (c[2] * (-sqrt(-4.0 * c[1] * c[3] + 4.0 *
                        c[1] * g_sol + c[2]^2) - c[2])^2) / (8.0 * c[1]^2) + (c[3] *
                        (-sqrt(-4.0 * c[1] * c[3] + 4.0 * c[1] * g_sol + c[2]^2) -
                        c[2])) / (2.0 * c[1]) - (g_sol * (-sqrt(-4.0 * c[1] * c[3] +
                        4.0 * c[1] * g_sol + c[2]^2) - c[2])) / (2.0 * c[1]))

                    if length(g_pump_vals[nw][a]) > 0
                        g_diff = maximum(abs.(g_pump_vals[nw][a] .- g_sol))
                    else
                        g_diff = Inf
                    end

                    if abs(g_nl_sol - g_nl_pred) > 1.0e-4 && g_diff > 1.0e-4
                        push!(g_pump_vals[nw][a], max(0.0, g_sol))
                        lhs = _calc_pump_gain_integrated_oa(g, z, g_sol, c)
                        con = JuMP.@build_constraint(lhs <= g_nl)
                        WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
                    end
                end
            end
        end
    end
end


function add_owf_user_cut_callback!(wm::WM.AbstractWaterModel)
    callback_function = get_owf_user_cut_callback(wm)
    WM._MOI.set(wm.model, WM._MOI.UserCutCallback(), callback_function)
end

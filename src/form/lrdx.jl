function _calc_pipe_flow_integrated_oa(q::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, q_hat::Float64, exponent::Float64)
    f = 1.0 / (1.0 + exponent) * q_hat^(1.0 + exponent) * z
    return f + q_hat^exponent * (q - q_hat * z)
end


function _calc_pipe_head_integrated_oa(dh::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, dh_hat::Float64, exponent::Float64)
    f = exponent / (1.0 + exponent) * dh_hat^(1.0 + 1.0 / exponent) * z
    return f + dh_hat^(1.0 / exponent) * (dh - dh_hat * z)
end


function _calc_pump_flow_integrated_oa(q::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, q_hat::Float64, coeffs::Array{Float64, 1})
    f = (1.0 / 3.0) * coeffs[1] * q_hat^3 + 0.5 * coeffs[2] * q_hat^2 + coeffs[3] * q_hat
    df = coeffs[1] * q_hat^2 + coeffs[2] * q_hat + coeffs[3]
    return f * z + df * (q - q_hat * z)
end


function _calc_pump_gain_integrated_oa(g::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, g_hat::Float64, c::Array{Float64, 1})
    f = (max(0.0, (4.0 * c[1] * (g_hat - c[3]) + c[2]^2))^(1.5) + c[2] *
        (6.0 * c[1] * (g_hat - c[3]) + c[2]^2)) / (12.0 * c[1]^2)
    df = (sqrt(max(0.0, 4.0 * c[1] * (g_hat - c[3]) + c[2]^2)) + c[2]) / (2.0 * c[1])
    return f * z + df * (g - g_hat * z)
end


function variable_pipe_flow_nonlinear(wm::AbstractLRDXModel; nw::Int=WM.nw_id_default, bounded::Bool=true, report::Bool=true)
    # Initialize variables associated with positive flows.
    qp_nl = WM.var(wm, nw)[:qp_nl_pipe] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound=0.0, base_name="$(nw)_qp_nl",
        start=WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "qp_nl_start", 1.0e-6))

    # Initialize variables associated with negative flows.
    qn_nl = WM.var(wm, nw)[:qn_nl_pipe] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound=0.0, base_name="$(nw)_qn_nl",
        start=WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "qn_nl_start", 1.0e-6))
end


function variable_pipe_head_difference_nonlinear(wm::AbstractLRDXModel; nw::Int=WM.nw_id_default, bounded::Bool=true, report::Bool=true)
    # Initialize variables associated with positive flows.
    dhp_nl = WM.var(wm, nw)[:dhp_nl_pipe] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound=0.0, base_name="$(nw)_dhp_nl",
        start=WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "dhp_nl_start", 1.0e-6))

    # Initialize variables associated with negative flows.
    dhn_nl = WM.var(wm, nw)[:dhn_nl_pipe] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound=0.0, base_name="$(nw)_dhn_nl",
        start=WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "dhn_nl_start", 1.0e-6))
end


function variable_pump_flow_nonlinear(wm::AbstractLRDXModel; nw::Int=WM.nw_id_default, bounded::Bool=true, report::Bool=true)
    # Initialize variables associated with positive flows.
    qp_nl = WM.var(wm, nw)[:qp_nl_pump] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pump)], lower_bound=0.0, base_name="$(nw)_qp_nl",
        start=WM.comp_start_value(WM.ref(wm, nw, :pump, a), "qp_nl_start", 1.0e-6))
end


function variable_pump_gain_nonlinear(wm::AbstractLRDXModel; nw::Int=WM.nw_id_default, bounded::Bool=true, report::Bool=true)
    # Initialize variables associated with positive gain.
    g_nl = WM.var(wm, nw)[:g_nl_pump] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pump)], lower_bound=0.0, base_name="$(nw)_g_nl",
        start=WM.comp_start_value(WM.ref(wm, nw, :pump, a), "g_nl_start", 1.0e-6))
end


function variable_tank_nonlinear(wm::AbstractLRDXModel; nw::Int=WM.nw_id_default, bounded::Bool=true, report::Bool=true)
    # Initialize variables associated with tank flow-head nonlinearities.
    qh_nl = WM.var(wm, nw)[:qh_nl_tank] = JuMP.@variable(
        wm.model, [i in WM.ids(wm, nw, :tank)], base_name="$(nw)_qh_nl",
        start=WM.comp_start_value(WM.ref(wm, nw, :tank, i), "qh_nl_start", 0.0))
end


function constraint_pipe_flow_nonlinear(
    wm::AbstractLRDXModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Get the number of breakpoints for the pipe.
    num_breakpoints = get(wm.ext, :pipe_breakpoints, 1)

    # Get the variable for flow directionality.
    y = WM.var(wm, n, :y_pipe, a)

    # Get variables for positive flow and head difference.
    qp, qp_nl = WM.var(wm, n, :qp_pipe, a), WM.var(wm, n, :qp_nl_pipe, a)
    qp_min_forward, qp_ub = max(0.0, q_min_forward), JuMP.upper_bound(qp)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(qp_min_forward, stop = JuMP.upper_bound(qp), length = num_breakpoints)
        # Add a linear outer approximation of the convex relaxation at `pt`.
        lhs = _calc_pipe_flow_integrated_oa(qp, y, pt, exponent)

        # Add outer-approximation of the integrated head loss constraint.
        c = JuMP.@constraint(wm.model, r * lhs <= qp_nl / L)

        # Append the :pipe_head_loss_integrated constraint array.
        append!(WM.con(wm, n, :pipe_flow_nonlinear)[a], [c])
    end

    # Get variables for negative flow and head difference.
    qn, qn_nl = WM.var(wm, n, :qn_pipe, a), WM.var(wm, n, :qn_nl_pipe, a)
    qn_min_forward, qn_ub = max(0.0, -q_max_reverse), JuMP.upper_bound(qn)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(qn_min_forward, stop = JuMP.upper_bound(qn), length = num_breakpoints)
        # Add a linear outer approximation of the convex relaxation at `pt`.
        lhs = _calc_pipe_flow_integrated_oa(qn, 1.0 - y, pt, exponent)

        # Add outer-approximation of the integrated head loss constraint.
        c = JuMP.@constraint(wm.model, r * lhs <= qn_nl / L)

        # Append the :pipe_head_loss_integrated constraint array.
        append!(WM.con(wm, n, :pipe_flow_nonlinear)[a], [c])
    end
end


function constraint_pipe_head_nonlinear(
    wm::AbstractLRDXModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Get the number of breakpoints for the pipe.
    num_breakpoints = get(wm.ext, :pipe_breakpoints, 1)

    # Get the variable for flow directionality.
    y = WM.var(wm, n, :y_pipe, a)
    L_r, r_r = L^(-1.0 / exponent), r^(-1.0 / exponent)

    # Get variables for positive flow and head difference.
    dhp, dhp_nl = WM.var(wm, n, :dhp_pipe, a), WM.var(wm, n, :dhp_nl_pipe, a)
    dhp_lb, dhp_ub = 0.0, JuMP.upper_bound(dhp)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(dhp_lb, stop = dhp_ub, length = num_breakpoints)
        # Add a linear outer approximation of the convex relaxation at `pt`.
        lhs = _calc_pipe_head_integrated_oa(dhp, y, pt, exponent)

        # Add outer-approximation of the integrated head loss constraint.
        c = JuMP.@constraint(wm.model, r_r * lhs <= dhp_nl / L_r)

        # Append the :pipe_head_loss_integrated constraint array.
        append!(WM.con(wm, n, :pipe_head_nonlinear)[a], [c])
    end

    # Get variables for positive flow and head difference.
    dhn, dhn_nl = WM.var(wm, n, :dhn_pipe, a), WM.var(wm, n, :dhn_nl_pipe, a)
    dhn_lb, dhn_ub = 0.0, JuMP.upper_bound(dhn)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(dhn_lb, stop = dhn_ub, length = num_breakpoints)
        # Add a linear outer approximation of the convex relaxation at `pt`.
        lhs = _calc_pipe_head_integrated_oa(dhn, 1.0 - y, pt, exponent)

        # Add outer-approximation of the integrated head loss constraint.
        c = JuMP.@constraint(wm.model, r_r * lhs <= dhn_nl / L_r)

        # Append the :pipe_head_loss_integrated constraint array.
        append!(WM.con(wm, n, :pipe_head_nonlinear)[a], [c])
    end
end


function constraint_on_off_pump_flow_nonlinear(
    wm::AbstractLRDXModel, n::Int, a::Int, node_fr::Int,
    node_to::Int, coeffs::Array{Float64, 1}, q_min_forward::Float64)
    # Get the number of breakpoints for the pump.
    num_breakpoints = get(wm.ext, :pump_breakpoints, 1)

    # Get the variable for pump status.
    z = WM.var(wm, n, :z_pump, a)

    # Get variables for positive flow and head difference.
    qp, qp_nl = WM.var(wm, n, :qp_pump, a), WM.var(wm, n, :qp_nl_pump, a)
    qp_min_forward, qp_ub = max(0.0, q_min_forward), JuMP.upper_bound(qp)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(qp_min_forward, stop = JuMP.upper_bound(qp), length = num_breakpoints)
        # Add a linear outer approximation of the convex relaxation at `pt`.
        rhs = _calc_pump_flow_integrated_oa(qp, z, pt, coeffs)

        # Add outer-approximation of the integrated head loss constraint.
        c = JuMP.@constraint(wm.model, qp_nl <= rhs)

        # Append the :pump_head_loss_integrated constraint array.
        append!(WM.con(wm, n, :on_off_pump_flow_nonlinear)[a], [c])
    end
end


function constraint_on_off_pump_gain_nonlinear(
    wm::AbstractLRDXModel, n::Int, a::Int, node_fr::Int,
    node_to::Int, coeffs::Array{Float64, 1}, q_min_forward::Float64)
    # Get the number of breakpoints for the pump.
    num_breakpoints = get(wm.ext, :pump_breakpoints, 1)

    # Get the variable for pump status.
    z = WM.var(wm, n, :z_pump, a)

    # Get variables for positive flow and head difference.
    g, g_nl = WM.var(wm, n, :g_pump, a), WM.var(wm, n, :g_nl_pump, a)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(0.0, stop = JuMP.upper_bound(g), length = num_breakpoints)
        # Add a linear outer approximation of the convex relaxation at `pt`.
        lhs = _calc_pump_gain_integrated_oa(g, z, pt, coeffs)

        # Add outer-approximation of the integrated head loss constraint.
        c = JuMP.@constraint(wm.model, lhs <= g_nl)

        # Append the :pump_head_loss_integrated constraint array.
        append!(WM.con(wm, n, :on_off_pump_gain_nonlinear)[a], [c])
    end
end


function constraint_tank_nonlinear(wm::AbstractLRDXModel, n::Int, i::Int, node_index::Int)
    q, h = WM.var(wm, n, :q_tank, i), WM.var(wm, n, :h, node_index)
    qh_nl_tank = WM.var(wm, n, :qh_nl_tank, i)

    q_lb, q_ub = JuMP.lower_bound(q), JuMP.upper_bound(q)
    h_lb, h_ub = JuMP.lower_bound(h), JuMP.upper_bound(h)

    c_1 = JuMP.@constraint(wm.model, qh_nl_tank >= q_lb * h + h_lb * q - q_lb * h_lb)
    c_2 = JuMP.@constraint(wm.model, qh_nl_tank >= q_ub * h + h_ub * q - q_ub * h_ub)
    c_3 = JuMP.@constraint(wm.model, qh_nl_tank <= q_lb * h + h_ub * q - q_lb * h_ub)
    c_4 = JuMP.@constraint(wm.model, qh_nl_tank <= q_ub * h + h_lb * q - q_ub * h_lb)

    append!(WM.con(wm, n, :tank_nonlinear)[i], [c_1, c_2, c_3, c_4])
end


""
function constraint_strong_duality(wm::AbstractLRDXModel)
    base_length = get(wm.data, "base_length", 1.0)
    base_time = get(wm.data, "base_time", 1.0)
    alpha = WM._get_alpha_min_1(wm) + 1.0
    
    pipe_type = wm.ref[:it][WM.wm_it_sym][:head_loss]
    viscosity = wm.ref[:it][WM.wm_it_sym][:viscosity]

    qp_pipe_nl = sum(sum(WM.var(wm, nw, :qp_nl_pipe)) for nw in WM.nw_ids(wm))
    qn_pipe_nl = sum(sum(WM.var(wm, nw, :qn_nl_pipe)) for nw in WM.nw_ids(wm))
    dhp_pipe_nl = sum(sum(WM.var(wm, nw, :dhp_nl_pipe)) for nw in WM.nw_ids(wm))
    dhn_pipe_nl = sum(sum(WM.var(wm, nw, :dhn_nl_pipe)) for nw in WM.nw_ids(wm))
    pipe_nl = qp_pipe_nl + qn_pipe_nl + dhp_pipe_nl + dhn_pipe_nl

    qp_pump_nl = sum(sum(WM.var(wm, nw, :qp_nl_pump)) for nw in WM.nw_ids(wm))
    g_pump_nl = sum(sum(WM.var(wm, nw, :g_nl_pump)) for nw in WM.nw_ids(wm))
    qh_nl_tank = sum(sum(WM.var(wm, nw, :qh_nl_tank)) for nw in WM.nw_ids(wm))
    pump_nl = qp_pump_nl - g_pump_nl

    reservoir_sum, demand_sum = JuMP.AffExpr(0.0), JuMP.AffExpr(0.0)

    for (n, network) in WM.nws(wm)
        for (i, res) in WM.ref(wm, n, :reservoir)
            head = WM.ref(wm, n, :node, res["node"])["head_nominal"]

            for comp in [:des_pipe, :pipe, :pump, :regulator, :short_pipe, :valve]
                qp = WM.var(wm, n, Symbol("qp_" * string(comp)))
                qn = WM.var(wm, n, Symbol("qn_" * string(comp)))

                for a in WM.ref(wm, n, Symbol(string(comp) * "_fr"), res["node"])
                    reservoir_sum += head * (qp[a] - qn[a])
                end

                for a in WM.ref(wm, n, Symbol(string(comp) * "_to"), res["node"])
                    reservoir_sum += head * (qn[a] - qp[a])
                end
            end
        end

        for (i, demand) in WM.ref(wm, n, :demand)
            demand_sum += demand["flow_nominal"] * WM.var(wm, n, :h, demand["node"])
        end
    end

    JuMP.@constraint(wm.model, pipe_nl - pump_nl + demand_sum - reservoir_sum + qh_nl_tank <= 0.0)
end
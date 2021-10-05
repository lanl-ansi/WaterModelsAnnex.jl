function _calc_pipe_flow_integrated_bound(q::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, q_lb::Float64, q_ub::Float64, exponent::Float64)
    f_lb, f_ub = q_lb^(1.0 + exponent), q_ub^(1.0 + exponent)
    return f_lb * z + (f_ub - f_lb) / (q_ub - q_lb) * (q - q_lb * z)
end


function _calc_pipe_flow_integrated_oa(q::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, q_hat::Float64, exponent::Float64)
    f = q_hat^(1.0 + exponent)
    df = (1.0 + exponent) * q_hat^exponent
    return f * z + df * (q - q_hat * z)
end


function _calc_pump_flow_integrated(q_hat::Float64, coeffs::Array{Float64, 1})
    return coeffs[1] * q_hat + coeffs[2] * q_hat^(1.0 + coeffs[3]) 
end


function _calc_pump_flow_integrated_oa(q::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, q_hat::Float64, coeffs::Array{Float64, 1})
    f = coeffs[1] * q_hat + coeffs[2] * q_hat^(1.0 + coeffs[3])
    df = coeffs[1] + (1.0 + coeffs[3]) * coeffs[2] * q_hat^(coeffs[3])
    return f * z + df * (q - q_hat * z)
end


function _calc_pump_flow_integrated_bound(q::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, q_lb::Float64, q_ub::Float64, coeffs::Array{Float64, 1})
    f_lb = coeffs[1] * q_lb + coeffs[2] * q_lb^(1.0 + coeffs[3])
    f_ub = coeffs[1] * q_ub + coeffs[2] * q_ub^(1.0 + coeffs[3])
    return f_lb * z + (f_ub - f_lb) / (q_ub - q_lb) * (q - q_lb * z)
end


function variable_pipe_flow_nonlinear(wm::Union{AbstractPWLRDXModel,AbstractLRDXModel}; nw::Int = WM.nw_id_default, bounded::Bool = true, report::Bool = true)
    # Initialize variables associated with positive flows.
    WM.var(wm, nw)[:qp_nl_pipe] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound = 0.0, base_name="$(nw)_qp_nl",
        start = WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "qp_nl_start", 1.0e-6))

    # Initialize variables associated with negative flows.
    WM.var(wm, nw)[:qn_nl_pipe] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound = 0.0, base_name = "$(nw)_qn_nl",
        start = WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "qn_nl_start", 1.0e-6))
end


function variable_pump_flow_nonlinear(wm::Union{AbstractPWLRDXModel,AbstractLRDXModel}; nw::Int = WM.nw_id_default, bounded::Bool = true, report::Bool = true)
    # Initialize variables associated with positive flows.
    WM.var(wm, nw)[:qp_nl_pump] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pump)], lower_bound = 0.0, base_name = "$(nw)_qp_nl",
        start = WM.comp_start_value(WM.ref(wm, nw, :pump, a), "qp_nl_start", 1.0e-6))
end


function variable_tank_nonlinear(wm::Union{AbstractPWLRDXModel,AbstractLRDXModel}; nw::Int=WM.nw_id_default, bounded::Bool = true, report::Bool = true)
    # Initialize variables associated with tank flow-head nonlinearities.
    WM.var(wm, nw)[:qh_nl_tank] = JuMP.@variable(
        wm.model, [i in WM.ids(wm, nw, :tank)], base_name = "$(nw)_qh_nl",
        start = WM.comp_start_value(WM.ref(wm, nw, :tank, i), "qh_nl_start", 0.0))
end


function constraint_pipe_flow_nonlinear(
    wm::AbstractLRDXModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Get the variable for flow directionality.
    y = WM.var(wm, n, :y_pipe, a)

    # Get variables for positive flow and nonlinear term.
    qp, qp_nl = WM.var(wm, n, :qp_pipe, a), WM.var(wm, n, :qp_nl_pipe, a)

    # Get the corresponding positive flow partitioning.
    partition_p = WM.get_pipe_flow_partition_positive(WM.ref(wm, n, :pipe, a))

    # Loop over consequential points (i.e., those that have nonzero head loss).
    for flow_value in filter(x -> x > 0.0, partition_p)
        # Add a linear outer approximation of the convex constraint at `pt`.
        lhs = _calc_pipe_flow_integrated_oa(qp, y, flow_value, exponent)

        # Add outer-approximation of the nonlinear flow constraint.
        c = JuMP.@constraint(wm.model, r * lhs <= qp_nl / L)

        # Append the :pipe_flow_nonlinear constraint array.
        append!(WM.con(wm, n, :pipe_flow_nonlinear)[a], [c])
    end

    # Get the corresponding min/max positive directed flows (when active).
    qp_min_forward, qp_max = max(0.0, q_min_forward), maximum(partition_p)

    if qp_min_forward != qp_max
        # Compute scaled versions of the nonlinear constraint overapproximations.
        f_1, f_2 = r * qp_min_forward^(1.0 + exponent), r * qp_max^(1.0 + exponent)
        f_slope = (f_2 - f_1) / (qp_max - qp_min_forward)
        f_lb_line = f_1 * y + f_slope * (qp - qp_min_forward * y)

        # Add upper-bounding lines of the head loss constraint.
        c = JuMP.@constraint(wm.model, qp_nl / L <= f_lb_line)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(WM.con(wm, n, :pipe_flow_nonlinear)[a], [c])
    end

    # Get variables for negative flow and nonlinear term.
    qn, qn_nl = WM.var(wm, n, :qn_pipe, a), WM.var(wm, n, :qn_nl_pipe, a)

    # Get the corresponding negative flow partitioning.
    partition_n = sort(-WM.get_pipe_flow_partition_negative(WM.ref(wm, n, :pipe, a)))

    # Loop over consequential points (i.e., those that have nonzero head loss).
    for flow_value in filter(x -> x > 0.0, partition_n)
        # Add a linear outer approximation of the convex relaxation at `flow_value`.
        lhs = _calc_pipe_flow_integrated_oa(qn, 1.0 - y, flow_value, exponent)

        # Add outer-approximation of the nonlinear flow constraint.
        c = JuMP.@constraint(wm.model, r * lhs <= qn_nl / L)

        # Append the :pipe_flow_nonlinear constraint array.
        append!(WM.con(wm, n, :pipe_flow_nonlinear)[a], [c])
    end

    # Get the corresponding maximum negative directed flow (when active).
    qn_min_forward, qn_max = max(0.0, -q_max_reverse), maximum(partition_n)

    if qn_min_forward != qn_max
        # Compute scaled versions of the head loss overapproximations.
        f_1, f_2 = r * qn_min_forward^(1.0 + exponent), r * qn_max^(1.0 + exponent)
        f_slope = (f_2 - f_1) / (qn_max - qn_min_forward)
        f_lb_line = f_slope * (qn - qn_min_forward * (1.0 - y)) + f_1 * (1.0 - y)

        # Add upper-bounding lines of the head loss constraint.
        c = JuMP.@constraint(wm.model, qn_nl / L <= f_lb_line)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(WM.con(wm, n, :pipe_flow_nonlinear)[a], [c])
    end
end


function constraint_on_off_pump_flow_nonlinear(
    wm::AbstractLRDXModel, n::Int, a::Int, node_fr::Int,
    node_to::Int, coeffs::Array{Float64, 1}, q_min_forward::Float64)
    # Get the variable for pump status.
    z = WM.var(wm, n, :z_pump, a)

    # Get variables for positive flow and head difference.
    qp, qp_nl = WM.var(wm, n, :qp_pump, a), WM.var(wm, n, :qp_nl_pump, a)
    partition = WM.ref(wm, n, :pump, a, "flow_partition")

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in partition
        # Add a linear outer approximation of the convex relaxation at `pt`.
        rhs = _calc_pump_flow_integrated_oa(qp, z, pt, coeffs)

        # Add outer-approximation of the integrated head loss constraint.
        c = JuMP.@constraint(wm.model, qp_nl <= rhs)

        # Append the :pump_head_loss_integrated constraint array.
        append!(WM.con(wm, n, :on_off_pump_flow_nonlinear)[a], [c])
    end

    # Get the corresponding min/max positive directed flows (when active).
    qp_min_forward, qp_max = max(0.0, q_min_forward), maximum(partition)

    if qp_min_forward != qp_max
        # Compute scaled versions of the head loss overapproximations.
        f_1 = _calc_pump_flow_integrated(qp_min_forward, coeffs)
        f_2 = _calc_pump_flow_integrated(qp_max, coeffs)        
        f_slope = (f_2 - f_1) / (qp_max - qp_min_forward)
        f_lb_line = f_1 * z + f_slope * (qp - qp_min_forward * z)

        # Add upper-bounding lines of the head loss constraint.
        c = JuMP.@constraint(wm.model, qp_nl <= f_lb_line)

        # Append the :on_off_des_pipe_head_loss constraint array.
        append!(WM.con(wm, n, :on_off_pump_flow_nonlinear)[a], [c])
    end
end


""
function constraint_tank_nonlinear(wm::Union{AbstractPWLRDXModel,AbstractLRDXModel}, n::Int, i::Int, node_index::Int)
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
function constraint_strong_duality(wm::Union{AbstractPWLRDXModel,AbstractLRDXModel}; nw::Int=WM.nw_id_default)
    qp_pipe_nl = sum(WM.var(wm, nw, :qp_nl_pipe))
    qn_pipe_nl = sum(WM.var(wm, nw, :qn_nl_pipe))
    qp_pump_nl = sum(WM.var(wm, nw, :qp_nl_pump))
    qh_nl_tank = sum(WM.var(wm, nw, :qh_nl_tank))

    reservoir_sum = JuMP.AffExpr(0.0)
    demand_sum = JuMP.AffExpr(0.0)

    for res in values(WM.ref(wm, nw, :reservoir))
        head = WM.ref(wm, nw, :node, res["node"])["head_nominal"]

        for comp in [:des_pipe, :pipe, :pump, :regulator, :short_pipe, :valve]
            qp = WM.var(wm, nw, Symbol("qp_" * string(comp)))
            qn = WM.var(wm, nw, Symbol("qn_" * string(comp)))

            for a in WM.ref(wm, nw, Symbol(string(comp) * "_fr"), res["node"])
                reservoir_sum += head * (qp[a] - qn[a])
            end

            for a in WM.ref(wm, nw, Symbol(string(comp) * "_to"), res["node"])
                reservoir_sum += head * (qn[a] - qp[a])
            end
        end
    end

    for demand in values(WM.ref(wm, nw, :demand))
        demand_sum += demand["flow_nominal"] * WM.var(wm, nw, :h, demand["node"])
    end

    linear_terms = demand_sum - reservoir_sum
    nonlinear_terms = qp_pipe_nl + qn_pipe_nl - qp_pump_nl - qh_nl_tank
    JuMP.@constraint(wm.model, linear_terms + nonlinear_terms <= 0.0)
end
function constraint_pipe_flow_nonlinear(
    wm::AbstractCDXModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Get the variable for flow directionality.
    y = WM.var(wm, n, :y_pipe, a)

    # Get variables for positive flow and nonlinear term.
    qp = WM.var(wm, n, :qp_pipe, a)
    qp_nl = WM.var(wm, n, :qp_nl_pipe, a)
    c_1 = JuMP.@NLconstraint(wm.model, r * qp^(exponent + 1.0) <= qp_nl / L)
    append!(WM.con(wm, n, :pipe_flow_nonlinear)[a], [c_1])

    # Get variables for negative flow and nonlinear term.
    qn = WM.var(wm, n, :qn_pipe, a)
    qn_nl = WM.var(wm, n, :qn_nl_pipe, a)
    c_2 = JuMP.@NLconstraint(wm.model, r * qn^(exponent + 1.0) <= qn_nl / L)
    append!(WM.con(wm, n, :pipe_flow_nonlinear)[a], [c_2])
end


function constraint_on_off_pump_flow_nonlinear(
    wm::AbstractCDXModel, n::Int, a::Int, node_fr::Int,
    node_to::Int, coeffs::Vector{Float64}, q_min_forward::Float64)
    qp = WM.var(wm, n, :qp_pump, a)
    qp_nl = WM.var(wm, n, :qp_nl_pump, a)
    c = JuMP.@NLconstraint(wm.model, qp_nl <= coeffs[1] * qp + coeffs[2] * qp^(1.0 + coeffs[3]))
    append!(WM.con(wm, n, :on_off_pump_flow_nonlinear)[a], [c])

    # # Get the variable for pump status.
    # z = WM.var(wm, n, :z_pump, a)

    # # Get variables for positive flow and head difference.
    # qp = WM.var(wm, n, :qp_pump, a)
    # qp_nl = WM.var(wm, n, :qp_nl_pump, a)
    # partition = WM.ref(wm, n, :pump, a, "flow_partition")

    # # Loop over breakpoints strictly between the lower and upper variable bounds.
    # for pt in partition
    #     # Add a linear outer approximation of the convex relaxation at `pt`.
    #     rhs = _calc_pump_flow_integrated_oa(qp, z, pt, coeffs)

    #     # Add outer-approximation of the integrated head loss constraint.
    #     c = JuMP.@constraint(wm.model, qp_nl <= rhs)

    #     # Append the :pump_head_loss_integrated constraint array.
    #     append!(WM.con(wm, n, :on_off_pump_flow_nonlinear)[a], [c])
    # end

    # # Get the corresponding min/max positive directed flows (when active).
    # qp_min_forward = max(0.0, q_min_forward)
    # qp_max = max(maximum(partition), JuMP.upper_bound(qp))

    # if qp_min_forward != qp_max
    #     # Add upper-bounding line for the nonlinear constraint.
    #     lhs = _calc_pump_flow_integrated_bound(qp, z, qp_min_forward, qp_max, coeffs)
    #     c = JuMP.@constraint(wm.model, lhs <= qp_nl)

    #     # Append the :on_off_des_pipe_head_loss constraint array.
    #     append!(WM.con(wm, n, :on_off_pump_flow_nonlinear)[a], [c])
    # end
end


""
function constraint_tank_nonlinear(wm::AbstractCDXModel, n::Int, i::Int, node_index::Int)
    q, h = WM.var(wm, n, :q_tank, i), WM.var(wm, n, :h, node_index)
    qh_nl_tank = WM.var(wm, n, :qh_nl_tank, i)
    JuMP.@constraint(wm.model, qh_nl_tank == q * h)
end
function constraint_pipe_flow_nonlinear(
    wm::AbstractPWLRDXModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Get the variable for flow directionality.
    y = WM.var(wm, n, :y_pipe, a)

    # Get variables for positive flow and nonlinear term.
    qp, qp_nl = WM.var(wm, n, :qp_pipe, a), WM.var(wm, n, :qp_nl_pipe, a)

    # Get the corresponding positive flow partitioning.
    partition_p = WM.get_pipe_flow_partition_positive(WM.ref(wm, n, :pipe, a))
    bp_range = 1:length(partition_p)

    # Loop over consequential points (i.e., those that have nonzero head loss).
    for flow_value in filter(x -> x > 0.0, partition_p)
        # Add a linear outer approximation of the convex relaxation at `pt`.
        lhs = _calc_pipe_flow_integrated_oa(qp, y, flow_value, exponent)

        # Add outer-approximation of the nonlinear flow constraint.
        c = JuMP.@constraint(wm.model, r * lhs <= qp_nl / L)

        # Append the :pipe_flow_nonlinear constraint array.
        append!(WM.con(wm, n, :pipe_flow_nonlinear)[a], [c])
    end

    # Add a constraint that upper-bounds the nonlinear flow variable.
    lambda_p = WM.var(wm, n, :lambda_p_pipe)
    f_p = r .* partition_p.^(1.0 + exponent)
    f_p_ub_expr = sum(f_p[k] * lambda_p[a, k] for k in bp_range)
    c = JuMP.@constraint(wm.model, qp_nl / L <= f_p_ub_expr)
    append!(WM.con(wm, n, :pipe_flow_nonlinear, a), [c])

    # Get variables for negative flow and nonlinear term.
    qn, qn_nl = WM.var(wm, n, :qn_pipe, a), WM.var(wm, n, :qn_nl_pipe, a)

    # Get the corresponding negative flow partitioning.
    partition_n = sort(-WM.get_pipe_flow_partition_negative(WM.ref(wm, n, :pipe, a)))
    bn_range = 1:length(partition_n)

    # Loop over consequential points (i.e., those that have nonzero head loss).
    for flow_value in filter(x -> x > 0.0, partition_n)
        # Add a linear outer approximation of the convex relaxation at `flow_value`.
        lhs = _calc_pipe_flow_integrated_oa(qn, 1.0 - y, flow_value, exponent)

        # Add outer-approximation of the nonlinear flow constraint.
        c = JuMP.@constraint(wm.model, r * lhs <= qn_nl / L)

        # Append the :pipe_flow_nonlinear constraint array.
        append!(WM.con(wm, n, :pipe_flow_nonlinear)[a], [c])
    end

    # Add a constraint that upper-bounds the head loss variable.
    lambda_n = WM.var(wm, n, :lambda_n_pipe)
    f_n = r .* partition_n.^(1.0 + exponent)
    f_n_ub_expr = sum(f_n[k] * lambda_n[a, k] for k in bn_range)
    c = JuMP.@constraint(wm.model, qn_nl / L <= f_n_ub_expr)
    append!(WM.con(wm, n, :pipe_flow_nonlinear, a), [c])
end


function constraint_on_off_pump_flow_nonlinear(
    wm::AbstractPWLRDXModel, n::Int, a::Int, node_fr::Int,
    node_to::Int, coeffs::Array{Float64, 1}, q_min_forward::Float64)
    # Get the variable for pump status.
    z = WM.var(wm, n, :z_pump, a)

    # Get variables for positive flow and nonlinear flow.
    qp, qp_nl = WM.var(wm, n, :qp_pump, a), WM.var(wm, n, :qp_nl_pump, a)
    partition = WM.ref(wm, n, :pump, a, "flow_partition")

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in partition
        # Add a linear outer approximation of the convex relaxation at `pt`.
        rhs = _calc_pump_flow_integrated_oa(qp, z, pt, coeffs)

        # Add outer-approximation of the integrated head gain constraint.
        c = JuMP.@constraint(wm.model, qp_nl <= rhs)

        # Append the :pump_head_loss_integrated constraint array.
        append!(WM.con(wm, n, :on_off_pump_flow_nonlinear)[a], [c])
    end

    # Add a constraint that lower-bounds the head gain variable.
    lambda = WM.var(wm, n, :lambda_pump)
    coeffs = WM.ref(wm, n, :pump, a, "head_curve_coefficients")
    f_all = _calc_pump_flow_integrated.(collect(partition), Ref(coeffs))
    f_lb_expr = sum(f_all[k] .* lambda[a, k] for k in 1:length(partition))
    c = JuMP.@constraint(wm.model, f_lb_expr <= qp_nl)
    append!(WM.con(wm, n, :on_off_pump_flow_nonlinear)[a], [c])
end
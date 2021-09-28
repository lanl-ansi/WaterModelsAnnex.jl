function _calc_pipe_flow_integrated_bound(q::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, q_lb::Float64, q_ub::Float64, exponent::Float64)
    f_lb = 1.0 / (1.0 + exponent) * q_lb^(1.0 + exponent)
    f_ub = 1.0 / (1.0 + exponent) * q_ub^(1.0 + exponent)
    return f_lb * z + (f_ub - f_lb) / (q_ub - q_lb) * (q - q_lb * z)
end


function _calc_pipe_flow_integrated_oa(q::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, q_hat::Float64, exponent::Float64)
    f = 1.0 / (1.0 + exponent) * q_hat^(1.0 + exponent) * z
    return f + q_hat^exponent * (q - q_hat * z)
end


function _calc_pipe_head_integrated_oa(dh::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, dh_hat::Float64, exponent::Float64)
    f = exponent / (1.0 + exponent) * dh_hat^(1.0 + 1.0 / exponent) * z
    return f + dh_hat^(1.0 / exponent) * (dh - dh_hat * z)
end


function _calc_pipe_head_integrated_bound(dh::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, dh_lb::Float64, dh_ub::Float64, exponent::Float64)
    f_lb = exponent / (1.0 + exponent) * dh_lb^(1.0 + 1.0 / exponent)
    f_ub = exponent / (1.0 + exponent) * dh_ub^(1.0 + 1.0 / exponent)
    return f_lb * z + (f_ub - f_lb) / (dh_ub - dh_lb) * (dh - dh_lb * z)
end


function _calc_pump_flow_integrated_oa(q::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, q_hat::Float64, coeffs::Array{Float64, 1})
    f = (1.0 / 3.0) * coeffs[1] * q_hat^3 + 0.5 * coeffs[2] * q_hat^2 + coeffs[3] * q_hat
    df = coeffs[1] * q_hat^2 + coeffs[2] * q_hat + coeffs[3]
    return f * z + df * (q - q_hat * z)
end


function _calc_pump_flow_integrated_bound(q::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, q_lb::Float64, q_ub::Float64, coeffs::Array{Float64, 1})
    f_lb = (1.0 / 3.0) * coeffs[1] * q_lb^3 + 0.5 * coeffs[2] * q_lb^2 + coeffs[3] * q_lb
    f_ub = (1.0 / 3.0) * coeffs[1] * q_ub^3 + 0.5 * coeffs[2] * q_ub^2 + coeffs[3] * q_ub
    return f_lb * z + (f_ub - f_lb) / (q_ub - q_lb) * (q - q_lb * z)
end


function _calc_pump_gain_integrated_oa(g::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, g_hat::Float64, c::Array{Float64, 1})
    f = ((-sqrt(-4.0 * c[1] * c[3] + 4.0 * c[1] * g_hat + c[2]^2) -
        c[2])^3 / (24.0 * c[1]^2) + (c[2] * (-sqrt(-4.0 * c[1] * c[3] + 4.0 *
        c[1] * g_hat + c[2]^2) - c[2])^2) / (8.0 * c[1]^2) + (c[3] *
        (-sqrt(-4.0 * c[1] * c[3] + 4.0 * c[1] * g_hat + c[2]^2) -
        c[2])) / (2.0 * c[1]) - (g_hat * (-sqrt(-4.0 * c[1] * c[3] +
        4.0 * c[1] * g_hat + c[2]^2) - c[2])) / (2.0 * c[1]))

    df = -(-sqrt(4.0 * c[1] * g_hat - 4.0 * c[1] * c[3] + c[2]^2) - c[2])^2 /
        (4.0 * c[1] * sqrt(4.0 * c[1] * g_hat - 4.0 * c[1] * c[3] + c[2]^2)) - 
        (c[2] * (-sqrt(4.0 * c[1] * g_hat - 4.0 * c[1] * c[3] + c[2]^2) - c[2])) /
        (2.0 * c[1] * sqrt(4.0 * c[1] * g_hat - 4.0 * c[1] * c[3] + c[2]^2)) -
        (-sqrt(4.0 * c[1] * g_hat - 4.0 * c[1] * c[3] + c[2]^2) - c[2]) /
        (2.0 * c[1]) + g_hat / sqrt(4.0 * c[1] * g_hat - 4.0 * c[1] * c[3] + c[2]^2) -
        c[3] / sqrt(4.0 * c[1] * g_hat - 4.0 * c[1] * c[3] + c[2]^2)

    return f * z + df * (g - g_hat * z)
end


function _calc_pump_gain_integrated_bound(g::JuMP.VariableRef, z::JuMP.VariableRef, g_lb::Float64, g_ub::Float64, c::Array{Float64, 1})
    f_lb = ((-sqrt(-4.0 * c[1] * c[3] + 4.0 * c[1] * g_lb + c[2]^2) -
        c[2])^3 / (24.0 * c[1]^2) + (c[2] * (-sqrt(-4.0 * c[1] * c[3] + 4.0 *
        c[1] * g_lb + c[2]^2) - c[2])^2) / (8.0 * c[1]^2) + (c[3] *
        (-sqrt(-4.0 * c[1] * c[3] + 4.0 * c[1] * g_lb + c[2]^2) -
        c[2])) / (2.0 * c[1]) - (g_lb * (-sqrt(-4.0 * c[1] * c[3] +
        4.0 * c[1] * g_lb + c[2]^2) - c[2])) / (2.0 * c[1]))

    f_ub = ((-sqrt(-4.0 * c[1] * c[3] + 4.0 * c[1] * g_ub + c[2]^2) -
        c[2])^3 / (24.0 * c[1]^2) + (c[2] * (-sqrt(-4.0 * c[1] * c[3] + 4.0 *
        c[1] * g_ub + c[2]^2) - c[2])^2) / (8.0 * c[1]^2) + (c[3] *
        (-sqrt(-4.0 * c[1] * c[3] + 4.0 * c[1] * g_ub + c[2]^2) -
        c[2])) / (2.0 * c[1]) - (g_ub * (-sqrt(-4.0 * c[1] * c[3] +
        4.0 * c[1] * g_ub + c[2]^2) - c[2])) / (2.0 * c[1]))

    return f_lb * z + (f_ub - f_lb) / (g_ub - g_lb) * (g - g_lb * z)
end


function variable_pipe_flow_nonlinear(wm::AbstractPWLRDXModel; nw::Int=WM.nw_id_default, bounded::Bool=true, report::Bool=true)
    # Initialize variables associated with positive flows.
    qp_nl = WM.var(wm, nw)[:qp_nl_pipe] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound=0.0, base_name="$(nw)_qp_nl",
        start=WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "qp_nl_start", 1.0e-6))

    # Initialize variables associated with negative flows.
    qn_nl = WM.var(wm, nw)[:qn_nl_pipe] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound=0.0, base_name="$(nw)_qn_nl",
        start=WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "qn_nl_start", 1.0e-6))
end


function variable_pipe_head_difference_nonlinear(wm::AbstractPWLRDXModel; nw::Int=WM.nw_id_default, bounded::Bool=true, report::Bool=true)
    # Initialize variables associated with positive flows.
    dhp_nl = WM.var(wm, nw)[:dhp_nl_pipe] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound=0.0, base_name="$(nw)_dhp_nl",
        start=WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "dhp_nl_start", 1.0e-6))

    # Initialize variables associated with negative flows.
    dhn_nl = WM.var(wm, nw)[:dhn_nl_pipe] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound=0.0, base_name="$(nw)_dhn_nl",
        start=WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "dhn_nl_start", 1.0e-6))
end


function variable_pump_flow_nonlinear(wm::AbstractPWLRDXModel; nw::Int=WM.nw_id_default, bounded::Bool=true, report::Bool=true)
    # Initialize variables associated with positive flows.
    qp_nl = WM.var(wm, nw)[:qp_nl_pump] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pump)], lower_bound=0.0, base_name="$(nw)_qp_nl",
        start=WM.comp_start_value(WM.ref(wm, nw, :pump, a), "qp_nl_start", 1.0e-6))
end


function variable_pump_gain_nonlinear(wm::AbstractPWLRDXModel; nw::Int=WM.nw_id_default, bounded::Bool=true, report::Bool=true)
    # Initialize variables associated with positive gain.
    g_nl = WM.var(wm, nw)[:g_nl_pump] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, :pump)], base_name="$(nw)_g_nl",
        start=WM.comp_start_value(WM.ref(wm, nw, :pump, a), "g_nl_start", 0.0))
end


function variable_tank_nonlinear(wm::AbstractPWLRDXModel; nw::Int=WM.nw_id_default, bounded::Bool=true, report::Bool=true)
    # Initialize variables associated with tank flow-head nonlinearities.
    qh_nl = WM.var(wm, nw)[:qh_nl_tank] = JuMP.@variable(
        wm.model, [i in WM.ids(wm, nw, :tank)], base_name="$(nw)_qh_nl",
        start=WM.comp_start_value(WM.ref(wm, nw, :tank, i), "qh_nl_start", 0.0))
end


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
    f_p = (r * 1.0 / (1.0 + exponent)) .* partition_p.^(1.0 + exponent)
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
    f_n = (r * 1.0 / (1.0 + exponent)) .* partition_n.^(1.0 + exponent)
    f_n_ub_expr = sum(f_n[k] * lambda_n[a, k] for k in bn_range)
    c = JuMP.@constraint(wm.model, qn_nl / L <= f_n_ub_expr)
    append!(WM.con(wm, n, :pipe_flow_nonlinear, a), [c])
end


function constraint_pipe_head_nonlinear(
    wm::AbstractPWLRDXModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Get object from the WaterModels reference dictionary.
    head_loss = wm.ref[:it][WM.wm_it_sym][:head_loss]
    viscosity = wm.ref[:it][WM.wm_it_sym][:viscosity]

    wm_data = WM.get_wm_data(wm.data)
    base_length = get(wm_data, "base_length", 1.0)
    base_mass = get(wm_data, "base_mass", 1.0)
    base_time = get(wm_data, "base_time", 1.0)

    # Get the variable for flow directionality.
    y = WM.var(wm, n, :y_pipe, a)
    L_r, r_r = L^(-1.0 / exponent), r^(-1.0 / exponent)

    # Get variables for positive head difference and nonlinear term.
    dhp, dhp_nl = WM.var(wm, n, :dhp_pipe, a), WM.var(wm, n, :dhp_nl_pipe, a)
    partition_p = WM.get_pipe_head_difference_partition_positive(
        WM.ref(wm, n, :pipe, a), head_loss, viscosity, base_length, base_mass, base_time)
    bp_range = 1:length(partition_p)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for dh_value in filter(x -> x > 0.0, partition_p)
        # Add a linear outer approximation of the convex relaxation at `dh_value`.
        lhs = _calc_pipe_head_integrated_oa(dhp, y, dh_value, exponent)

        # Add outer-approximation of the integrated head loss constraint.
        c = JuMP.@constraint(wm.model, r_r * lhs <= dhp_nl / L_r)

        # Append the :pipe_head_nonlinear constraint array.
        append!(WM.con(wm, n, :pipe_head_nonlinear)[a], [c])
    end

    # Add a constraint that upper-bounds the nonlinear flow variable.
    lambda_p = WM.var(wm, n, :lambda_p_pipe)
    f_p = (r_r * exponent / (1.0 + exponent)) .* partition_p.^(1.0 + 1.0 / exponent)
    f_p_ub_expr = sum(f_p[k] * lambda_p[a, k] for k in bp_range)
    c = JuMP.@constraint(wm.model, dhp_nl / L_r <= f_p_ub_expr)
    append!(WM.con(wm, n, :pipe_head_nonlinear)[a], [c])

    # Get variables for positive flow and head difference.
    dhn, dhn_nl = WM.var(wm, n, :dhn_pipe, a), WM.var(wm, n, :dhn_nl_pipe, a)
    partition_n = sort(WM.get_pipe_head_difference_partition_negative(
        WM.ref(wm, n, :pipe, a), head_loss, viscosity, base_length, base_mass, base_time))
    bn_range = 1:length(partition_n)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for dh_value in filter(x -> x > 0.0, partition_n)
        # Add a linear outer approximation of the convex relaxation at `dh_value`.
        lhs = _calc_pipe_head_integrated_oa(dhn, 1.0 - y, dh_value, exponent)

        # Add outer-approximation of the integrated head loss constraint.
        c = JuMP.@constraint(wm.model, r_r * lhs <= dhn_nl / L_r)

        # Append the :pipe_head_nonlinear constraint array.
        append!(WM.con(wm, n, :pipe_head_nonlinear)[a], [c])
    end

    # Add a constraint that upper-bounds the nonlinear flow variable.
    lambda_n = WM.var(wm, n, :lambda_n_pipe)
    f_n = (r_r * exponent / (1.0 + exponent)) .* partition_n.^(1.0 + 1.0 / exponent)
    f_n_ub_expr = sum(f_n[k] * lambda_n[a, k] for k in bn_range)
    c = JuMP.@constraint(wm.model, dhn_nl / L_r <= f_n_ub_expr)
    append!(WM.con(wm, n, :pipe_head_nonlinear)[a], [c])
end


function constraint_on_off_pump_flow_nonlinear(
    wm::AbstractPWLRDXModel, n::Int, a::Int, node_fr::Int,
    node_to::Int, coeffs::Array{Float64, 1}, q_min_forward::Float64)
    # Get object from the WaterModels reference dictionary.
    pump = WM.ref(wm, n, :pump, a)

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

    # Add a constraint that lower-bounds the head gain variable.
    lambda, g = WM.var(wm, n, :lambda_pump), WM.var(wm, n, :g_pump, a)
    head_curve_function = WM._calc_head_curve_function(pump)
    f_all = head_curve_function.(collect(partition))
    gain_lb_expr = sum(f_all[k] .* lambda[a, k] for k in 1:length(partition))
    c = JuMP.@constraint(wm.model, gain_lb_expr <= g)
    append!(WM.con(wm, n, :on_off_pump_flow_nonlinear)[a], [c])
end


function constraint_on_off_pump_gain_nonlinear(
    wm::AbstractPWLRDXModel, n::Int, a::Int, node_fr::Int,
    node_to::Int, coeffs::Array{Float64, 1}, q_min_forward::Float64)
    # Get object from the WaterModels reference dictionary.
    pump = WM.ref(wm, n, :pump, a)
 
    # Get the variable for pump status.
    z = WM.var(wm, n, :z_pump, a)

    # Get variables for positive flow and head difference.
    g = WM.var(wm, n, :g_pump, a)
    g_nl = WM.var(wm, n, :g_nl_pump, a)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in WM.get_pump_head_gain_partition(pump)
        # Add a linear outer approximation of the convex relaxation at `pt`.
        lhs = _calc_pump_gain_integrated_oa(g, z, pt, coeffs)

        # Add outer-approximation of the integrated head gain constraint.
        c = JuMP.@constraint(wm.model, lhs <= g_nl)

        # Append the :on_off_pump_gain_nonlinear constraint array.
        append!(WM.con(wm, n, :on_off_pump_gain_nonlinear)[a], [c])
    end

    lambda, f_all = WM.var(wm, n, :lambda_pump), Vector{Float64}([])
    breakpoints = WM.get_pump_head_gain_partition(pump)

    for g_hat in breakpoints
        push!(f_all, ((-sqrt(-4.0 * coeffs[1] * coeffs[3] + 4.0 * coeffs[1] * g_hat + coeffs[2]^2) -
            coeffs[2])^3 / (24.0 * coeffs[1]^2) + (coeffs[2] * (-sqrt(-4.0 * coeffs[1] * coeffs[3] + 4.0 *
            coeffs[1] * g_hat + coeffs[2]^2) - coeffs[2])^2) / (8.0 * coeffs[1]^2) + (coeffs[3] *
            (-sqrt(-4.0 * coeffs[1] * coeffs[3] + 4.0 * coeffs[1] * g_hat + coeffs[2]^2) -
            coeffs[2])) / (2.0 * coeffs[1]) - (g_hat * (-sqrt(-4.0 * coeffs[1] * coeffs[3] +
            4.0 * coeffs[1] * g_hat + coeffs[2]^2) - coeffs[2])) / (2.0 * coeffs[1])))
    end

    gain_lb_expr = sum(f_all[k] .* lambda[a, k] for k in 1:length(breakpoints))
    c = JuMP.@constraint(wm.model, g_nl <= gain_lb_expr)
    append!(WM.con(wm, n, :on_off_pump_flow_nonlinear)[a], [c])
end


function constraint_tank_nonlinear(wm::AbstractPWLRDXModel, n::Int, i::Int, node_index::Int)
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
function constraint_strong_duality(wm::AbstractPWLRDXModel; nw::Int=WM.nw_id_default)
    qp_pipe_nl = sum(WM.var(wm, nw, :qp_nl_pipe))
    qn_pipe_nl = sum(WM.var(wm, nw, :qn_nl_pipe))
    dhp_pipe_nl = sum(WM.var(wm, nw, :dhp_nl_pipe))
    dhn_pipe_nl = sum(WM.var(wm, nw, :dhn_nl_pipe))
    pipe_nl = qp_pipe_nl + qn_pipe_nl + dhp_pipe_nl + dhn_pipe_nl

    qp_pump_nl = sum(WM.var(wm, nw, :qp_nl_pump))
    g_pump_nl = sum(WM.var(wm, nw, :g_nl_pump))
    qh_nl_tank = sum(WM.var(wm, nw, :qh_nl_tank))
    pump_nl = qp_pump_nl - g_pump_nl

    reservoir_sum, demand_sum = JuMP.AffExpr(0.0), JuMP.AffExpr(0.0)

    for (i, res) in WM.ref(wm, nw, :reservoir)
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

    for (i, demand) in WM.ref(wm, nw, :demand)
        demand_sum += demand["flow_nominal"] * WM.var(wm, nw, :h, demand["node"])
    end

    linear_terms = demand_sum - reservoir_sum
    nonlinear_terms = pipe_nl - pump_nl - qh_nl_tank

    JuMP.@constraint(wm.model, linear_terms + nonlinear_terms <= 0.0)
end
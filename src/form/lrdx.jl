function _calc_pipe_flow_integrated_oa(q::JuMP.VariableRef, z::Union{JuMP.VariableRef, JuMP.GenericAffExpr}, q_hat::Float64, exponent::Float64)
    f = 1.0 / (1.0 + exponent) * q_hat^(1.0 + exponent) * z
    return f + q_hat^exponent * (q - q_hat * z)
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


function constraint_pipe_flow_nonlinear(
    wm::AbstractLRDXModel, n::Int, a::Int, node_fr::Int, node_to::Int, exponent::Float64,
    L::Float64, r::Float64, q_max_reverse::Float64, q_min_forward::Float64)
    # Get the number of breakpoints for the pipe.
    num_breakpoints = get(wm.ext, :pipe_breakpoints, 1)

    # Get the variable for flow directionality.
    y = WM.var(wm, n, :y_pipe, a)

    # Get variables for positive flow and head difference.
    qp, qp_nl = WM.var(wm, n, :qp_pipe, a), WM.var(wm, n, :qp_nl_pipe, a)
    dhp_lb, dhp_ub = 0.0, JuMP.upper_bound(dhp)
    qp_min_forward, qp_ub = max(0.0, q_min_forward), JuMP.upper_bound(qp)

    # Loop over breakpoints strictly between the lower and upper variable bounds.
    for pt in range(qp_min_forward, stop = JuMP.upper_bound(qp), length = num_breakpoints+2)[2:end-1]
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
    for pt in range(qn_min_forward, stop = JuMP.upper_bound(qn), length = num_breakpoints+2)[2:end-1]
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
    println(num_breakpoints)

    # # Get the variable for flow directionality.
    # y = WM.var(wm, n, :y_pipe, a)

    # # Get variables for positive flow and head difference.
    # qp, qp_nl = WM.var(wm, n, :qp_pipe, a), WM.var(wm, n, :qp_nl_pipe, a)
    # qp_min_forward, qp_ub = max(0.0, q_min_forward), JuMP.upper_bound(qp)

    # # Loop over breakpoints strictly between the lower and upper variable bounds.
    # for pt in range(qp_min_forward, stop = JuMP.upper_bound(qp), length = num_breakpoints+2)[2:end-1]
    #     # Add a linear outer approximation of the convex relaxation at `pt`.
    #     lhs = _calc_pipe_flow_integrated_oa(qp, y, pt, exponent)

    #     # Add outer-approximation of the integrated head loss constraint.
    #     c = JuMP.@constraint(wm.model, r * lhs <= qp_nl / L)

    #     # Append the :pipe_head_loss_integrated constraint array.
    #     append!(WM.con(wm, n, :pipe_head_nonlinear)[a], [c])
    # end

    # # Get variables for negative flow and head difference.
    # qn, qn_nl = WM.var(wm, n, :qn_pipe, a), WM.var(wm, n, :qn_nl_pipe, a)
    # qn_min_forward, qn_ub = max(0.0, -q_max_reverse), JuMP.upper_bound(qn)

    # # Loop over breakpoints strictly between the lower and upper variable bounds.
    # for pt in range(qn_min_forward, stop = JuMP.upper_bound(qn), length = num_breakpoints+2)[2:end-1]
    #     # Add a linear outer approximation of the convex relaxation at `pt`.
    #     lhs = _calc_pipe_flow_integrated_oa(qn, 1.0 - y, pt, exponent)

    #     # Add outer-approximation of the integrated head loss constraint.
    #     c = JuMP.@constraint(wm.model, r * lhs <= qn_nl / L)

    #     # Append the :pipe_head_loss_integrated constraint array.
    #     append!(WM.con(wm, n, :pipe_head_nonlinear)[a], [c])
    # end
end


"""
    objective_wf(wm::AbstractLRDXModel)
"""
function constraint_strong_duality(wm::AbstractLRDXModel)
    base_length = get(wm.data, "base_length", 1.0)
    base_time = get(wm.data, "base_time", 1.0)
    alpha = WM._get_alpha_min_1(wm) + 1.0
    
    pipe_type = wm.ref[:it][WM.wm_it_sym][:head_loss]
    viscosity = wm.ref[:it][WM.wm_it_sym][:viscosity]


    qp_pipe_nl = sum(sum(WM.var(wm, nw, :qp_nl_pipe)) for nw in WM.nw_ids(wm))
    qn_pipe_nl = sum(sum(WM.var(wm, nw, :qn_nl_pipe)) for nw in WM.nw_ids(wm))

    println(qp_pipe_nl)

    # q_pipe_nl_sum = sum()

    # f_1, f_2 = Array{Any, 1}([0.0]), Array{Any, 1}([0.0])
    # f_3, f_4 = Array{Any, 1}([0.0]), Array{Any, 1}([0.0])

    # for (nw, network) in WM.nws(wm)
    #     # Initialize variables associated with positive flows.
    #     qp_pipe_nl = 

    #     # Initialize variables associated with negative flows.
    #     qn_pipe_nl = WM.var(wm, nw)[:qn_pipe_nl] = JuMP.@variable(
    #         wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound=0.0,base_name="$(nw)_qn",
    #         start=WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "qn_nl_start", WM._FLOW_MIN))
    # end



    # # f_5, f_6 = Array{Any, 1}([0.0]), Array{Any, 1}([0.0])

    # for (nw, network) in WM.nws(wm)
    #     # Initialize variables associated with positive flows.
    #     qp_pipe_nl = WM.var(wm, nw)[:qp_pipe_nl] = JuMP.@variable(
    #         wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound=0.0, base_name="$(nw)_qp_nl",
    #         start=WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "qp_nl_start", WM._FLOW_MIN))

    #     # Initialize variables associated with negative flows.
    #     qn_pipe_nl = WM.var(wm, nw)[:qn_pipe_nl] = JuMP.@variable(
    #         wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound=0.0,base_name="$(nw)_qn",
    #         start=WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "qn_nl_start", WM._FLOW_MIN))
    
    #     # Initialize variables associated with positive flows.
    #     dhp_pipe_nl = WM.var(wm, nw)[:dhp_pipe_nl] = JuMP.@variable(
    #         wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound=0.0, base_name="$(nw)_dhp_nl",
    #         start=WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "dhp_nl_start", 1.0e-6))
    
    #     # Initialize variables associated with negative flows.
    #     dhn_pipe_nl = WM.var(wm, nw)[:dhn_pipe_nl] = JuMP.@variable(
    #         wm.model, [a in WM.ids(wm, nw, :pipe)], lower_bound=0.0,base_name="$(nw)_dhn",
    #         start=WM.comp_start_value(WM.ref(wm, nw, :pipe, a), "dhn_nl_start", 1.0e-6))

    #     # Initialize variables associated with positive flows.
    #     qp_pump_nl = WM.var(wm, nw)[:qp_pump_nl] = JuMP.@variable(
    #         wm.model, [a in WM.ids(wm, nw, :pump)], lower_bound=0.0, base_name="$(nw)_qp_nl",
    #         start=WM.comp_start_value(WM.ref(wm, nw, :pump, a), "qp_nl_start", WM._FLOW_MIN))

    #     # Get head variables.
    #     h = WM.var(wm, n, :h)

    #     # Get pipe flow variables.
    #     qp_pipe = WM.var(wm, n, :qp_pipe)
    #     qn_pipe = WM.var(wm, n, :qn_pipe)

    #     # Get pipe head difference variables.
    #     dhp_pipe = WM.var(wm, n, :dhp_pipe)
    #     dhn_pipe = WM.var(wm, n, :dhn_pipe)

    #     # Get design pipe flow variables.
    #     qp_des_pipe = WM.var(wm, n, :qp_des_pipe)
    #     qn_des_pipe = WM.var(wm, n, :qn_des_pipe)

    #     # Get design pipe head difference variables.
    #     dhp_des_pipe = WM.var(wm, n, :dhp_des_pipe)
    #     dhn_des_pipe = WM.var(wm, n, :dhn_des_pipe)

    #     # Get valve flow variables.
    #     qp_valve = WM.var(wm, n, :qp_valve)
    #     qn_valve = WM.var(wm, n, :qn_valve)

    #     # Get pump flow and head difference variables.
    #     q_tank = WM.var(wm, n, :q_tank)
    #     qp_pump = WM.var(wm, n, :qp_pump)
    #     g_pump = WM.var(wm, n, :g_pump)
    #     z_pump = WM.var(wm, n, :z_pump)
        
    #     for (a, pipe) in WM.ref(wm, n, :pipe)
    #         L_x_r = pipe["length"] * WM._calc_pipe_resistance(pipe, pipe_type, viscosity, base_length, base_time)
    #         push!(f_1, JuMP.@NLexpression(wm.model, L_x_r * head_loss(qp_pipe[a])))
    #         push!(f_1, JuMP.@NLexpression(wm.model, L_x_r * head_loss(qn_pipe[a])))
    #         push!(f_3, JuMP.@NLexpression(wm.model, L_x_r^(-1.0 / alpha) * head_loss_dh(dhp_pipe[a])))
    #         push!(f_3, JuMP.@NLexpression(wm.model, L_x_r^(-1.0 / alpha) * head_loss_dh(dhn_pipe[a])))
    #     end

    # #     for (a, des_pipe) in WM.ref(wm, n, :des_pipe)
    # #         L_x_r = des_pipe["length"] * WM._calc_pipe_resistance(des_pipe, pipe_type, viscosity, base_length, base_time)
    # #         push!(f_1, JuMP.@NLexpression(wm.model, L_x_r * head_loss(qp_des_pipe[a])))
    # #         push!(f_1, JuMP.@NLexpression(wm.model, L_x_r * head_loss(qn_des_pipe[a])))
    # #         push!(f_3, JuMP.@NLexpression(wm.model, L_x_r^(-1.0 / alpha) * head_loss_dh(dhp_des_pipe[a])))
    # #         push!(f_3, JuMP.@NLexpression(wm.model, L_x_r^(-1.0 / alpha) * head_loss_dh(dhn_des_pipe[a])))
    # #     end

    # #     for (a, pump) in WM.ref(wm, n, :pump)
    # #         @assert pump["head_curve_form"] in [WM.PUMP_QUADRATIC, WM.PUMP_BEST_EFFICIENCY_POINT, WM.PUMP_LINEAR_POWER]
    # #         c = WM._calc_head_curve_coefficients(pump)
            
    # #         push!(f_5, JuMP.@NLexpression(wm.model, c[1] / 3.0 * qp_pump[a]^3 +
    # #             0.5 * c[2] * qp_pump[a]^2 + c[3] * qp_pump[a]))

    # #         push!(f_6, JuMP.@NLexpression(wm.model, z_pump[a] * (
    # #             (-sqrt(-4.0 * c[1] * c[3] + 4.0 * c[1] * g_pump[a] + c[2]^2) -
    # #             c[2])^3 / (24.0 * c[1]^2) + (c[2] * (-sqrt(-4.0 * c[1] * c[3] + 4.0 *
    # #             c[1] * g_pump[a] + c[2]^2) - c[2])^2) / (8.0 * c[1]^2) + (c[3] *
    # #             (-sqrt(-4.0 *c[1] * c[3] + 4.0 * c[1] * g_pump[a] + c[2]^2) -
    # #             c[2])) / (2.0 * c[1]) - (g_pump[a] * (-sqrt(-4.0 * c[1] * c[3] +
    # #             4.0 * c[1] * g_pump[a] + c[2]^2) - c[2])) / (2.0 * c[1]))))
    # #     end

    # #     for (i, res) in WM.ref(wm, n, :reservoir)
    # #         head = WM.ref(wm, n, :node, res["node"])["head_nominal"]

    # #         for comp in [:des_pipe, :pipe, :pump, :regulator, :short_pipe, :valve]
    # #             qp = WM.var(wm, n, Symbol("qp_" * string(comp)))
    # #             qn = WM.var(wm, n, Symbol("qn_" * string(comp)))

    # #             for a in WM.ref(wm, n, Symbol(string(comp) * "_fr"), res["node"])
    # #                 push!(f_2, JuMP.@NLexpression(wm.model, head * (qp[a] - qn[a])))
    # #             end

    # #             for a in WM.ref(wm, n, Symbol(string(comp) * "_to"), res["node"])
    # #                 push!(f_2, JuMP.@NLexpression(wm.model, head * (qn[a] - qp[a])))
    # #             end
    # #         end
    # #     end

    # #     for (i, demand) in WM.ref(wm, n, :demand)
    # #         push!(f_4, JuMP.@NLexpression(wm.model, demand["flow_nominal"] * h[demand["node"]]))
    # #     end

    # #     for (i, tank) in WM.ref(wm, n, :tank)
    # #         # TODO: How should we convexify this?
    # #         head = tank["init_level"] + WM.ref(wm, n, :node, tank["node"])["elevation"]
    # #         push!(f_4, JuMP.@NLexpression(wm.model, -q_tank[i] * head))
    # #     end
    # end

    # JuMP.@NLconstraint(wm.model,
    #     sum(f_1[k] for k in 1:length(f_1)) -
    #     sum(f_2[k] for k in 1:length(f_2)) +
    #     sum(f_3[k] for k in 1:length(f_3)) +
    #     sum(f_4[k] for k in 1:length(f_4)) -
    #     sum(f_5[k] for k in 1:length(f_5)) +
    #     sum(f_6[k] for k in 1:length(f_6)) <= 0.0)
end
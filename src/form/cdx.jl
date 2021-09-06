function variable_flow_nonlinear(wm::AbstractCDXModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    
end


"Create flow-related variables common to all directed flow models for node-connecting components."
function variable_flow(wm::AbstractCDXModel; nw::Int=nw_id_default, bounded::Bool=true, report::Bool=true)
    for name in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
        # Create directed flow (`qp` and `qn`) variables for each component.
        _variable_component_flow(wm, name; nw=nw, bounded=bounded, report=report)

        # Create directed flow binary direction variables (`y`) for each component.
        _variable_component_direction(wm, name; nw=nw, report=report)
    end

    for name in ["des_pipe", "pipe"]
        # Create directed head difference (`dhp` and `dhn`) variables for each component.
        _variable_component_head_difference(wm, name; nw=nw, bounded=bounded, report=report)
    end
end


"""
    objective_wf(wm::AbstractCDXModel)
"""
function constraint_strong_duality(wm::AbstractCDXModel)
    base_length = get(wm.data, "base_length", 1.0)
    base_time = get(wm.data, "base_time", 1.0)
    alpha = WM._get_alpha_min_1(wm) + 1.0
    
    pipe_type = wm.ref[:it][WM.wm_it_sym][:head_loss]
    viscosity = wm.ref[:it][WM.wm_it_sym][:viscosity]

    f_1, f_2 = Array{Any, 1}([0.0]), Array{Any, 1}([0.0])
    f_3, f_4 = Array{Any, 1}([0.0]), Array{Any, 1}([0.0])
    f_5, f_6 = Array{Any, 1}([0.0]), Array{Any, 1}([0.0])

    for (n, network) in WM.nws(wm)
        # Get head variables.
        h = WM.var(wm, n, :h)

        # Get pipe flow variables.
        qp_pipe = WM.var(wm, n, :qp_pipe)
        qn_pipe = WM.var(wm, n, :qn_pipe)

        # Get pipe head difference variables.
        dhp_pipe = WM.var(wm, n, :dhp_pipe)
        dhn_pipe = WM.var(wm, n, :dhn_pipe)

        # Get design pipe flow variables.
        qp_des_pipe = WM.var(wm, n, :qp_des_pipe)
        qn_des_pipe = WM.var(wm, n, :qn_des_pipe)

        # Get design pipe head difference variables.
        dhp_des_pipe = WM.var(wm, n, :dhp_des_pipe)
        dhn_des_pipe = WM.var(wm, n, :dhn_des_pipe)

        # Get valve flow variables.
        qp_valve = WM.var(wm, n, :qp_valve)
        qn_valve = WM.var(wm, n, :qn_valve)

        # Get pump flow and head difference variables.
        q_tank = WM.var(wm, n, :q_tank)
        qp_pump = WM.var(wm, n, :qp_pump)
        g_pump = WM.var(wm, n, :g_pump)
        z_pump = WM.var(wm, n, :z_pump)
        
        for (a, pipe) in WM.ref(wm, n, :pipe)
            L_x_r = pipe["length"] * WM._calc_pipe_resistance(pipe, pipe_type, viscosity, base_length, base_time)
            push!(f_1, JuMP.@NLexpression(wm.model, L_x_r * head_loss(qp_pipe[a])))
            push!(f_1, JuMP.@NLexpression(wm.model, L_x_r * head_loss(qn_pipe[a])))
            push!(f_3, JuMP.@NLexpression(wm.model, L_x_r^(-1.0 / alpha) * head_loss_dh(dhp_pipe[a])))
            push!(f_3, JuMP.@NLexpression(wm.model, L_x_r^(-1.0 / alpha) * head_loss_dh(dhn_pipe[a])))
        end

        for (a, des_pipe) in WM.ref(wm, n, :des_pipe)
            L_x_r = des_pipe["length"] * WM._calc_pipe_resistance(des_pipe, pipe_type, viscosity, base_length, base_time)
            push!(f_1, JuMP.@NLexpression(wm.model, L_x_r * head_loss(qp_des_pipe[a])))
            push!(f_1, JuMP.@NLexpression(wm.model, L_x_r * head_loss(qn_des_pipe[a])))
            push!(f_3, JuMP.@NLexpression(wm.model, L_x_r^(-1.0 / alpha) * head_loss_dh(dhp_des_pipe[a])))
            push!(f_3, JuMP.@NLexpression(wm.model, L_x_r^(-1.0 / alpha) * head_loss_dh(dhn_des_pipe[a])))
        end

        for (a, pump) in WM.ref(wm, n, :pump)
            @assert pump["pump_type"] in [WM.PUMP_QUADRATIC, WM.PUMP_BEST_EFFICIENCY_POINT, WM.PUMP_LINEAR_POWER]
            c = WM._calc_head_curve_coefficients(pump)
            
            push!(f_5, JuMP.@NLexpression(wm.model, c[1] / 3.0 * qp_pump[a]^3 +
                0.5 * c[2] * qp_pump[a]^2 + c[3] * qp_pump[a]))

            push!(f_6, JuMP.@NLexpression(wm.model, z_pump[a] * (
                (-sqrt(-4.0 * c[1] * c[3] + 4.0 * c[1] * g_pump[a] + c[2]^2) -
                c[2])^3 / (24.0 * c[1]^2) + (c[2] * (-sqrt(-4.0 * c[1] * c[3] + 4.0 *
                c[1] * g_pump[a] + c[2]^2) - c[2])^2) / (8.0 * c[1]^2) + (c[3] *
                (-sqrt(-4.0 *c[1] * c[3] + 4.0 * c[1] * g_pump[a] + c[2]^2) -
                c[2])) / (2.0 * c[1]) - (g_pump[a] * (-sqrt(-4.0 * c[1] * c[3] +
                4.0 * c[1] * g_pump[a] + c[2]^2) - c[2])) / (2.0 * c[1]))))
        end

        for (i, res) in WM.ref(wm, n, :reservoir)
            head = WM.ref(wm, n, :node, res["node"])["head_nominal"]

            for comp in [:des_pipe, :pipe, :pump, :regulator, :short_pipe, :valve]
                qp = WM.var(wm, n, Symbol("qp_" * string(comp)))
                qn = WM.var(wm, n, Symbol("qn_" * string(comp)))

                for a in WM.ref(wm, n, Symbol(string(comp) * "_fr"), res["node"])
                    push!(f_2, JuMP.@NLexpression(wm.model, head * (qp[a] - qn[a])))
                end

                for a in WM.ref(wm, n, Symbol(string(comp) * "_to"), res["node"])
                    push!(f_2, JuMP.@NLexpression(wm.model, head * (qn[a] - qp[a])))
                end
            end
        end

        for (i, demand) in WM.ref(wm, n, :demand)
            push!(f_4, JuMP.@NLexpression(wm.model, demand["flow_nominal"] * h[demand["node"]]))
        end

        for (i, tank) in WM.ref(wm, n, :tank)
            # TODO: How should we convexify this?
            head = tank["init_level"] + WM.ref(wm, n, :node, tank["node"])["elevation"]
            push!(f_4, JuMP.@NLexpression(wm.model, -q_tank[i] * head))
        end
    end

    JuMP.@NLconstraint(wm.model,
        sum(f_1[k] for k in 1:length(f_1)) -
        sum(f_2[k] for k in 1:length(f_2)) +
        sum(f_3[k] for k in 1:length(f_3)) +
        sum(f_4[k] for k in 1:length(f_4)) -
        sum(f_5[k] for k in 1:length(f_5)) +
        sum(f_6[k] for k in 1:length(f_6)) <= 0.0)
end
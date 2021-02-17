# Define common CD (convex directed) implementations of water distribution network constraints.


"Create flow variables that are common to all directed flow models for a component."
function WM._variable_component_flow(
    wm::AbstractCDModel, component_name::String; nw::Int=wm.cnw,
    bounded::Bool=true, report::Bool=true)
    # Store the corresponding component symbol.
    comp_sym = Symbol(component_name)

    # Initialize variables associated with positive flows.
    qp = WM.var(wm, nw)[Symbol("qp_" * component_name)] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, comp_sym)], lower_bound=0.0, base_name="$(nw)_qp",
        start=WM.comp_start_value(WM.ref(wm, nw, comp_sym, a), "qp_start", WM._FLOW_MIN))

    # Initialize variables associated with negative flows.
    qn = WM.var(wm, nw)[Symbol("qn_" * component_name)] = JuMP.@variable(
        wm.model, [a in WM.ids(wm, nw, comp_sym)], lower_bound=0.0, base_name="$(nw)_qn",
        start=WM.comp_start_value(WM.ref(wm, nw, comp_sym, a), "qn_start", WM._FLOW_MIN))

    if bounded # Bound flow-related variables if desired.
        for (a, comp) in WM.ref(wm, nw, comp_sym)
            qp_max = max(0.0, comp["flow_max"])
            JuMP.set_upper_bound(qp[a], qp_max)
            qp_start = WM.comp_start_value(comp, "qp_start", 0.5 * qp_max)
            JuMP.set_start_value(qp[a], qp_start)

            qn_max = max(0.0, -comp["flow_min"])
            JuMP.set_upper_bound(qn[a], qn_max)
            qn_start = WM.comp_start_value(comp, "qn_start", 0.0)
            JuMP.set_start_value(qn[a], qn_start)
        end
    end

    # Report positive directed flow values as part of the solution.
    report && WM.sol_component_value(wm, nw, comp_sym, :qp, WM.ids(wm, nw, comp_sym), qp)

    # Report negative directed flow values as part of the solution.
    report && WM.sol_component_value(wm, nw, comp_sym, :qn, WM.ids(wm, nw, comp_sym), qn)

    # Create expressions capturing the relationships among q, qp, and qn.
    q = WM.var(wm, nw)[Symbol("q_" * component_name)] = JuMP.@expression(
        wm.model, [a in WM.ids(wm, nw, comp_sym)], qp[a] - qn[a])

    # Report flow expression values as part of the solution.
    report && WM.sol_component_value(wm, nw, comp_sym, :q, WM.ids(wm, nw, comp_sym), q)
end


"Create flow-related variables common to all directed flow models for node-connecting components."
function WM.variable_flow(wm::AbstractCDModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    for name in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
        # Create directed flow (`qp` and `qn`) variables for each component.
        WM._variable_component_flow(wm, name; nw=nw, bounded=bounded, report=report)
    end

    for name in ["des_pipe", "pipe"]
        # Create directed head difference (`dhp` and `dhn`) variables for each component.
        WM._variable_component_head_difference(wm, name; nw=nw, bounded=bounded, report=report)
    end
end


function WM.constraint_pipe_head(wm::AbstractCDModel, n::Int, a::Int, node_fr::Int, node_to::Int)
    # Get common flow variables and associated data.
    dhp, dhn = WM.var(wm, n, :dhp_pipe, a), WM.var(wm, n, :dhn_pipe, a)

    # Get head variables for from and to nodes.
    h_i, h_j = WM.var(wm, n, :h, node_fr), WM.var(wm, n, :h, node_to)

    # For pipes, the differences must satisfy lower and upper bounds.
    c = JuMP.@constraint(wm.model, dhp - dhn == h_i - h_j)

    # Append the :pipe_head constraint array.
    append!(WM.con(wm, n, :pipe_head, a), [c])
end


function WM.constraint_on_off_des_pipe_flow(wm::AbstractCDModel, n::Int, a::Int, q_max_reverse::Float64, q_min_forward::Float64)
    # Get des_pipe status variable.
    qp, qn = WM.var(wm, n, :qp_des_pipe, a), WM.var(wm, n, :qn_des_pipe, a)
    z = WM.var(wm, n, :z_des_pipe, a)

    # If the des_pipe is inactive, flow must be zero.
    qp_ub, qn_ub = JuMP.upper_bound(qp), JuMP.upper_bound(qn)
    c_1 = JuMP.@constraint(wm.model, qp <= qp_ub * z)
    c_2 = JuMP.@constraint(wm.model, qn <= qn_ub * z)

    # Append the constraint array.
    append!(WM.con(wm, n, :on_off_des_pipe_flow, a), [c_1, c_2])
end


function WM.constraint_on_off_des_pipe_head(wm::AbstractCDModel, n::Int, a::Int, node_fr::Int, node_to::Int)
    # Get head difference variables for the des_pipe.
    dhp, dhn = WM.var(wm, n, :dhp_des_pipe, a), WM.var(wm, n, :dhn_des_pipe, a)

    # Get des_pipe status variable.
    z = WM.var(wm, n, :z_des_pipe, a)

    # If the des_pipe is off, decouple the head difference relationship.
    dhp_ub, dhn_ub = JuMP.upper_bound(dhp), JuMP.upper_bound(dhn)
    c_1 = JuMP.@constraint(wm.model, dhp <= dhp_ub * z)
    c_2 = JuMP.@constraint(wm.model, dhn <= dhn_ub * z)

    # Append the constraint array.
    append!(WM.con(wm, n, :on_off_des_pipe_head, a), [c_1, c_2])
end


function WM.constraint_on_off_pump_flow(wm::AbstractCDModel, n::Int, a::Int, q_min_forward::Float64)
    # Get pump status variable.
    qp, z = WM.var(wm, n, :qp_pump, a), WM.var(wm, n, :z_pump, a)

    # If the pump is inactive, flow must be zero.
    qp_lb, qp_ub = q_min_forward, JuMP.upper_bound(qp)
    c_1 = JuMP.@constraint(wm.model, qp >= qp_lb * z)
    c_2 = JuMP.@constraint(wm.model, qp <= qp_ub * z)

    # Append the constraint array.
    append!(WM.con(wm, n, :on_off_pump_flow, a), [c_1, c_2])
end


function WM.constraint_on_off_regulator_flow(wm::AbstractCDModel, n::Int, a::Int, q_min_forward::Float64)
    # Get regulator flow, status, and direction variables.
    qp, z = WM.var(wm, n, :qp_regulator, a), WM.var(wm, n, :z_regulator, a)

    # If the regulator is closed, flow must be zero.
    qp_lb, qp_ub = max(JuMP.lower_bound(qp), q_min_forward), JuMP.upper_bound(qp)
    c_1 = JuMP.@constraint(wm.model, qp >= qp_lb * z)
    c_2 = JuMP.@constraint(wm.model, qp <= qp_ub * z)

    # Append the :on_off_regulator_flow constraint array.
    append!(WM.con(wm, n, :on_off_regulator_flow, a), [c_1, c_2])
end


function WM.constraint_short_pipe_flow(wm::AbstractCDModel, n::Int, a::Int, q_max_reverse::Float64, q_min_forward::Float64)
    # Nothing here for an AbstractCDModel.
end


function WM.constraint_on_off_valve_flow(wm::AbstractCDModel, n::Int, a::Int, q_max_reverse::Float64, q_min_forward::Float64)
    # Get valve flow, direction, and status variables.
    qp, qn = WM.var(wm, n, :qp_valve, a), WM.var(wm, n, :qn_valve, a)
    z = WM.var(wm, n, :z_valve, a)

    # The valve flow is constrained by direction and status.
    qp_ub, qn_ub = JuMP.upper_bound(qp), JuMP.upper_bound(qn)
    c_1 = JuMP.@constraint(wm.model, qp <= qp_ub * z)
    c_2 = JuMP.@constraint(wm.model, qn <= qn_ub * z)

    # Append the constraint array.
    append!(WM.con(wm, n, :on_off_valve_flow, a), [c_1, c_2])
end


"""
    objective_wf(wm::AbstractCDModel)
"""
function constraint_strong_duality(wm::AbstractCDModel)
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

        # Get pump flow and head difference variables.
        q_tank = WM.var(wm, n, :q_tank)
        qp_pump = WM.var(wm, n, :qp_pump)
        g_pump = WM.var(wm, n, :g_pump)
        
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
            @assert pump["head_curve_form"] in [WM.QUADRATIC, WM.BEST_EFFICIENCY_POINT, WM.LINEAR_POWER]
            c = WM._calc_head_curve_coefficients(pump)
            
            push!(f_5, JuMP.@NLexpression(wm.model, c[1] / 3.0 * qp_pump[a]^3 +
                0.5 * c[2] * qp_pump[a]^2 + c[3] * qp_pump[a]))

            push!(f_6, JuMP.@NLexpression(wm.model,
                (-sqrt(-4.0 * c[1] * c[3] + 4.0 * c[1] * g_pump[a] + c[2]^2) -
                c[2])^3 / (24.0 * c[1]^2) + (c[2] * (-sqrt(-4.0 * c[1] * c[3] + 4.0 *
                c[1] * g_pump[a] + c[2]^2) - c[2])^2) / (8.0 * c[1]^2) + (c[3] *
                (-sqrt(-4.0 *c[1] * c[3] + 4.0 * c[1] * g_pump[a] + c[2]^2) -
                c[2])) / (2.0 * c[1]) - (g_pump[a] * (-sqrt(-4.0 * c[1] * c[3] +
                4.0 * c[1] * g_pump[a] + c[2]^2) - c[2])) / (2.0 * c[1])))
        end

        for (i, res) in WM.ref(wm, n, :reservoir)
            head = WM.ref(wm, n, :node, res["node"])["head_nominal"]

            for a in WM.ref(wm, n, :pipe_fr, res["node"])
                push!(f_2, JuMP.@NLexpression(wm.model, head * (qp_pipe[a] - qn_pipe[a])))
            end

            for a in WM.ref(wm, n, :pipe_to, res["node"])
                push!(f_2, JuMP.@NLexpression(wm.model, head * (qn_pipe[a] - qp_pipe[a])))
            end

            for a in WM.ref(wm, n, :pump_fr, res["node"])
                push!(f_2, JuMP.@NLexpression(wm.model, head * qp_pump[a]))
            end

            for a in WM.ref(wm, n, :pump_to, res["node"])
                push!(f_2, JuMP.@NLexpression(wm.model, -head * qp_pump[a]))
            end

            for a in WM.ref(wm, n, :des_pipe_fr, res["node"])
                push!(f_2, JuMP.@NLexpression(wm.model, head * (qp_des_pipe[a] - qn_des_pipe[a])))
            end

            for a in WM.ref(wm, n, :des_pipe_to, res["node"])
                push!(f_2, JuMP.@NLexpression(wm.model, head * (qn_des_pipe[a] - qp_des_pipe[a])))
            end
        end

        for (i, demand) in WM.ref(wm, n, :demand)
            push!(f_4, JuMP.@NLexpression(wm.model, demand["flow_nominal"] * h[demand["node"]]))
        end

        for (i, tank) in WM.ref(wm, n, :tank)
            # TODO: How should we convexify this?
            push!(f_4, JuMP.@NLexpression(wm.model, -q_tank[i] * h[tank["node"]]))
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
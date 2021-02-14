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
function variable_flow(wm::AbstractCDModel; nw::Int=wm.cnw, bounded::Bool=true, report::Bool=true)
    for name in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
        # Create directed flow (`qp` and `qn`) variables for each component.
        WM._variable_component_flow(wm, name; nw=nw, bounded=bounded, report=report)
    end

    for name in ["des_pipe", "pipe"]
        # Create directed head difference (`dhp` and `dhn`) variables for each component.
        WM._variable_component_head_difference(wm, name; nw=nw, bounded=bounded, report=report)
    end
end

"""
    objective_wf(wm::AbstractCDModel)
"""
function WM.objective_wf(wm::AbstractCDModel)
    n = wm.cnw

    qp_pipe, qn_pipe = WM.var(wm, n, :qp_pipe), WM.var(wm, n, :qn_pipe)

    pipe_type = wm.ref[:it][WM.wm_it_sym][:head_loss]
    viscosity = wm.ref[:it][WM.wm_it_sym][:viscosity]

    L_x_r = Dict{Int, Any}(a => pipe["length"] *
        WM._calc_pipe_resistance(pipe, pipe_type, viscosity)
        for (a, pipe) in WM.ref(wm, n, :pipe))

    f_1 = JuMP.@NLexpression(wm.model,
        sum(L_x_r[a] * 0.35063113604 *
        (qp_pipe[a]^2.852 + qn_pipe[a]^2.852)
        for (a, pipe) in WM.ref(wm, n, :pipe)))

    f_2 = JuMP.@NLexpression(wm.model,
        sum(WM.ref(wm, n, :node, reservoir["node"])["head_nominal"]
        * sum(qp_pipe[a] - qn_pipe[a] for (a) in
        WM.ref(wm, n, :pipe_fr, reservoir["node"]))
        for (i, reservoir) in WM.ref(wm, n, :reservoir)))

    JuMP.@NLobjective(wm.model, MOI.MIN_SENSE, f_1 - f_2)
end
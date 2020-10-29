function compute_tau(wm::GenericWaterModel, q::Float64,
    mu_e::Float64, h::Dict{Int, Float64}, mu_v::Dict{Int, Float64},
    lambda_e_p::Float64, lambda_e_n::Float64, zeta::Float64,
    a::Int, optimizer::JuMP.OptimizerFactory, res::Float64, 
    q_lb::Float64, q_ub::Float64, n::Int=wm.cnw)
    L = WM.ref(wm, n, :links, a)["length"]
    i = WM.ref(wm, n, :links, a)["node_fr"]
    j = WM.ref(wm, n, :links, a)["node_to"]

    model = JuMP.Model(with_optimizer(optimizer))
    q_e = JuMP.@variable(model, lower_bound=q_lb, upper_bound=q_ub, base_name="q")
    alpha = WM.get_alpha_min_1(wm)
    JuMP.register(model, :head_loss, 1, WM.f_alpha(alpha, convex=false),
        WM.df_alpha(alpha, convex=false), WM.d2f_alpha(alpha, convex=false))

    obj = JuMP.@NLobjective(model, MOI.MIN_SENSE, (zeta * (q_e - q) + (1 - zeta) * mu_e) *
        (L * res * head_loss(q_e)) - zeta * (h[i] - h[j]) * q_e - (1 - zeta) *
        (mu_v[i] - mu_v[j] - lambda_e_p + lambda_e_n) * q_e)

    JuMP.optimize!(model)
    return JuMP.objective_value(model)
end

function compute_zeta(wm::GenericWaterModel, q::Dict{Int, Float64},
    mu_e::Dict{Int, Float64}, h::Dict{Int, Float64}, mu_v::Dict{Int, Float64},
    lambda_e_p::Dict{Int, Float64}, lambda_e_n::Dict{Int, Float64},
    optimizer::JuMP.OptimizerFactory, n::Int=wm.cnw)
    model = JuMP.Model(with_optimizer(optimizer))
    zeta = JuMP.@variable(model, base_name="zeta", lower_bound=0.0, upper_bound=1.0)

    for a in keys(q)
        if mu_e[a] * q[a] > 0.0
            JuMP.@constraint(model, zeta * abs(q[a]) >= (1 - zeta) * abs(mu_e[a]))
        elseif mu_e[a] * q[a] < 0.0
            i = WM.ref(wm, n, :links, a)["node_fr"]
            j = WM.ref(wm, n, :links, a)["node_to"]
            JuMP.@constraint(model, (1 - zeta) * abs(mu_v[i] - mu_v[j] - lambda_e_p[a] + lambda_e_n[a]) <= zeta * abs(mu_v[i] - mu_v[j]))
        else
            JuMP.@constraint(model, (1 - zeta) * mu_e[a] == 0.0)
        end
    end

    JuMP.@objective(model, MOI.MIN_SENSE, zeta)
    JuMP.optimize!(model)
    return JuMP.objective_value(model)
end

function add_humpola_cut!(wm::GenericWaterModel, optimizer::JuMP.OptimizerFactory, n::Int=wm.cnw)
    q_sol, h_sol, resistance_indices = get_cnlp_solution(wm, optimizer)
    h_ids = sort(collect(keys(h_sol)))
    q_ids = sort(collect(keys(q_sol)))
    alpha = WM.ref(wm, n, :alpha)

    model = JuMP.Model(with_optimizer(optimizer))
    q = JuMP.@variable(model, [a in q_ids], base_name="q")
    h = JuMP.@variable(model, [i in h_ids], base_name="h")
    qr = JuMP.@variable(model, [i in ids(wm, n, :reservoirs)], base_name="qr")
    Delta_e = JuMP.@variable(model, [a in q_ids], lower_bound=0.0, base_name="Delta_e")
    Delta_v = JuMP.@variable(model, [i in h_ids], lower_bound=0.0, base_name="Delta_v")

    loss = Dict{Int, JuMP.ConstraintRef}()
    conservation = Dict{Int, JuMP.ConstraintRef}()
    head_ub = Dict{Int, JuMP.ConstraintRef}()
    head_lb = Dict{Int, JuMP.ConstraintRef}()
    flow_lb = Dict{Int, JuMP.ConstraintRef}()
    flow_ub = Dict{Int, JuMP.ConstraintRef}()

    alpha = WM.get_alpha_min_1(wm)
    JuMP.register(model, :head_loss, 1, WM.f_alpha(alpha, convex=false),
        WM.df_alpha(alpha, convex=false), WM.d2f_alpha(alpha, convex=false))

    for a in q_ids
        r = resistance_indices[a]
        i = WM.ref(wm, n, :links, a)["node_fr"]
        j = WM.ref(wm, n, :links, a)["node_to"]
        L = WM.ref(wm, n, :links, a)["length"]
        res = WM.ref(wm, n, :resistance, a)[r]
        loss[a] = JuMP.@NLconstraint(model, L*res * head_loss(q[a]) == h[i] - h[j])

        q_ub = JuMP.upper_bound(WM.var(wm, n, :qp_ne, a)[r])
        q_lb = -JuMP.upper_bound(WM.var(wm, n, :qn_ne, a)[r])
        flow_ub[a] = JuMP.@constraint(model, q[a] - Delta_e[a] <= q_ub)
        flow_lb[a] = JuMP.@constraint(model, q[a] + Delta_e[a] >= q_lb)
    end

    for i in h_ids
        arcs_fr = ref(wm, n, :node_arcs_fr, i)
        arcs_to = ref(wm, n, :node_arcs_to, i)
        res = ref(wm, n, :node_reservoirs, i)
        junc = ref(wm, n, :node_junctions, i)
        demands = Dict(k => ref(wm, n, :junctions, k, "demand") for k in junc)

        c = JuMP.@constraint(model,
            sum( q[l] for (l, f, t) in arcs_fr) +
            sum(-q[l] for (l, f, t) in arcs_to) ==
            sum( qr[rid] for rid in res) +
            sum(-demand for (jid, demand) in demands))

        conservation[i] = c

        h_ub = JuMP.upper_bound(WM.var(wm, n, :h, i))
        h_lb = JuMP.lower_bound(WM.var(wm, n, :h, i))
        head_ub[i] = JuMP.@constraint(model, h[i] - Delta_v[i] <= h_ub)
        head_lb[i] = JuMP.@constraint(model, h[i] + Delta_v[i] >= h_lb)
    end

    JuMP.@objective(model, MOI.MIN_SENSE, sum(Delta_v) + sum(Delta_e))
    JuMP.optimize!(model)

    lambda_v_p = Dict{Int, Float64}(i => max(0, -JuMP.dual(c)) for (i, c) in head_ub)
    lambda_v_n = Dict{Int, Float64}(i => max(0, JuMP.dual(c)) for (i, c) in head_lb)
    lambda_v = Dict{Int, Float64}(i => max(0, JuMP.dual(JuMP.LowerBoundRef(Delta_v[i]))) for i in h_ids)
    lambda_e_p = Dict{Int, Float64}(a => max(0, -JuMP.dual(c)) for (a, c) in flow_ub)
    lambda_e_n = Dict{Int, Float64}(a => max(0, JuMP.dual(c)) for (a, c) in flow_lb)
    lambda_e = Dict{Int, Float64}(a => max(0, JuMP.dual(JuMP.LowerBoundRef(Delta_e[a]))) for a in q_ids)

    mu_v = Dict{Int, Float64}(i => JuMP.dual(c) for (i, c) in conservation)
    mu_e = Dict{Int, Float64}(a => -JuMP.dual(c) for (a, c) in loss)
    h = Dict{Int, Float64}(i => JuMP.value(h[i]) for i in h_ids)
    q = Dict{Int, Float64}(a => JuMP.value(q[a]) for a in q_ids)
    qr = Dict{Int, Float64}(i => JuMP.value(qr[i]) for i in ids(wm, n, :reservoirs))
    Delta_v = Dict{Int, Float64}(i => JuMP.value(Delta_v[i]) for i in h_ids)
    Delta_e = Dict{Int, Float64}(a => JuMP.value(Delta_e[a]) for a in q_ids)
    zeta = compute_zeta(wm, q, mu_e, h, mu_v, lambda_e_p, lambda_e_n, optimizer, n)

    lhs = JuMP.AffExpr(0.0)
    rhs = JuMP.AffExpr(0.0)

    rhs -= zeta * sum(h[i] * -junc["demand"] for (i, junc) in ref(wm, n, :junctions))
    rhs -= zeta * sum(h[i] * qr[i] for (i, res) in ref(wm, n, :reservoirs))
    rhs -= (1 - zeta) * sum(mu_v[i] * -junc["demand"] for (i, junc) in ref(wm, n, :junctions))
    rhs -= (1 - zeta) * sum(mu_v[i] * qr[i] for (i, res) in ref(wm, n, :reservoirs))
    rhs += (1 - zeta) * sum(lambda_v_p[i] * JuMP.upper_bound(WM.var(wm, n, :h, i)) for i in h_ids)
    rhs -= (1 - zeta) * sum(lambda_v_n[i] * JuMP.lower_bound(WM.var(wm, n, :h, i)) for i in h_ids)

    for a in q_ids
        for (r_id, r) in enumerate(ref(wm, n, :resistance, a))
            i = WM.ref(wm, n, :links, a)["node_fr"]
            j = WM.ref(wm, n, :links, a)["node_to"]
            res = WM.ref(wm, n, :resistance, a)[r_id]

            h_i_ub = JuMP.upper_bound(WM.var(wm, n, :h, i))
            h_i_lb = JuMP.lower_bound(WM.var(wm, n, :h, i))
            h_j_ub = JuMP.upper_bound(WM.var(wm, n, :h, j))
            h_j_lb = JuMP.lower_bound(WM.var(wm, n, :h, j))

            q_ub = JuMP.upper_bound(WM.var(wm, n, :qp_ne, a)[r_id])
            q_lb = -JuMP.upper_bound(WM.var(wm, n, :qn_ne, a)[r_id])

            mu_e_sol = r_id == resistance_indices[a] ? mu_e[a] : 0.0
            lambda_e_p_sol = r_id == resistance_indices[a] ? lambda_e_p[a] : 0.0
            lambda_e_n_sol = r_id == resistance_indices[a] ? lambda_e_n[a] : 0.0

            tau = compute_tau(wm, q[a], mu_e_sol, h, mu_v, lambda_e_p_sol, lambda_e_n_sol, zeta, a, optimizer, res, q_lb, q_ub, n)
            lhs += tau * var(wm, n, :x_res, a)[r_id]

            if r_id == resistance_indices[a]
                rhs += (1 - zeta) * var(wm, n, :x_res, a)[r_id] * (lambda_e_p_sol * q_ub)
                rhs -= (1 - zeta) * var(wm, n, :x_res, a)[r_id] * (lambda_e_n_sol * q_lb)
            end
        end
    end

    JuMP.@constraint(wm.model, lhs <= rhs)
end

function compute_tau(wm::GenericWaterModel, q::Dict{Int, Float64},
    mu_e::Dict{Int, Float64}, h::Dict{Int, Float64}, mu_v::Dict{Int, Float64},
    lambda_e_p::Dict{Int, Float64}, lambda_e_n::Dict{Int, Float64}, zeta::Float64,
    a::Int, optimizer::JuMP.OptimizerFactory, res::Float64, n::Int=wm.cnw)
    L = WMs.ref(wm, n, :links, a)["length"]

    i = WMs.ref(wm, n, :links, a)["node_fr"]
    j = WMs.ref(wm, n, :links, a)["node_to"]

    model = JuMP.Model(with_optimizer(optimizer))
    q_e = JuMP.@variable(model, base_name="q")
    alpha = WMs.get_alpha_min_1(wm)
    JuMP.register(model, :head_loss, 1, WMs.f_alpha(alpha, convex=false),
        WMs.df_alpha(alpha, convex=false), WMs.d2f_alpha(alpha, convex=false))

    obj = JuMP.@NLobjective(model, MOI.MIN_SENSE, (zeta * (q_e - q[a]) + (1 - zeta) * mu_e[a]) *
        (L * res * head_loss(q_e)) - zeta * (h[i] - h[j]) * q_e - (1 - zeta) *
        (mu_v[i] - mu_v[j] - lambda_e_p[a] + lambda_e_n[a]) * q_e)

    JuMP.optimize!(model)
    return JuMP.objective_value(model)
end

function compute_zeta(wm::GenericWaterModel, q::Dict{Int, Float64},
    mu_e::Dict{Int, Float64}, h::Dict{Int, Float64}, mu_v::Dict{Int, Float64},
    lambda_e_p::Dict{Int, Float64}, lambda_e_n::Dict{Int, Float64}, n::Int=wm.cnw)

    converged = false
    zeta_m1 = 0.0
    zeta = 0.0

    while !converged
        for a in keys(q)
            zeta_m1 = zeta

            i = WMs.ref(wm, n, :links, a)["node_fr"]
            j = WMs.ref(wm, n, :links, a)["node_to"]

            if mu_e[a] * q[a] > 0.0
                if !(zeta * abs(q[a]) > (1 - zeta) * abs(mu_e[a]))
                    zeta += 1.0e-5
                    break
                end
            end

            if mu_e[a] * q[a] < 0.0
                if !((1 - zeta) * abs(mu_v[i] - mu_v[j] - lambda_e_p[a] + lambda_e_n[a]) < zeta * abs(mu_v[i] - mu_v[j]))
                    zeta += 1.0e-5
                    break
                end
            end
        end

        if isapprox(zeta, zeta_m1)
            converged = true
        end
    end

    return zeta
end

function add_humpola_cut!(wm::GenericWaterModel, optimizer::JuMP.OptimizerFactory, n::Int=wm.cnw)
    q_sol, h_sol, resistance_indices = get_cnlp_solution(wm, optimizer)
    h_ids = sort(collect(keys(h_sol)))
    q_ids = sort(collect(keys(q_sol)))
    alpha = WMs.ref(wm, n, :alpha)

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

    alpha = WMs.get_alpha_min_1(wm)
    JuMP.register(model, :head_loss, 1, WMs.f_alpha(alpha, convex=false),
        WMs.df_alpha(alpha, convex=false), WMs.d2f_alpha(alpha, convex=false))

    for a in q_ids
        r = resistance_indices[a]
        i = WMs.ref(wm, n, :links, a)["node_fr"]
        j = WMs.ref(wm, n, :links, a)["node_to"]
        L = WMs.ref(wm, n, :links, a)["length"]
        res = WMs.ref(wm, n, :resistance, a)[r]
        loss[a] = JuMP.@NLconstraint(model, L*res * head_loss(q[a]) == h[i] - h[j])

        q_ub = JuMP.upper_bound(WMs.var(wm, n, :qp_ne, a)[r])
        q_lb = -JuMP.upper_bound(WMs.var(wm, n, :qn_ne, a)[r])
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
            sum(-q[l] for (l, f, t) in arcs_fr) +
            sum( q[l] for (l, f, t) in arcs_to) ==
            sum(-qr[rid] for rid in res) +
            sum(demand for (jid, demand) in demands))

        conservation[i] = c

        h_ub = JuMP.upper_bound(WMs.var(wm, n, :h, i))
        h_lb = JuMP.lower_bound(WMs.var(wm, n, :h, i))
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

    #lambda_v_p = Dict{Int, Float64}(i => JuMP.dual(c) for (i, c) in head_ub)
    #lambda_v_n = Dict{Int, Float64}(i => JuMP.dual(c) for (i, c) in head_lb)
    #lambda_v = Dict{Int, Float64}(i => JuMP.dual(JuMP.LowerBoundRef(Delta_v[i])) for i in h_ids)
    #lambda_e_p = Dict{Int, Float64}(a => JuMP.dual(c) for (a, c) in flow_ub)
    #lambda_e_n = Dict{Int, Float64}(a => JuMP.dual(c) for (a, c) in flow_lb)
    #lambda_e = Dict{Int, Float64}(a => JuMP.dual(JuMP.LowerBoundRef(Delta_e[a])) for a in q_ids)

    mu_v = Dict{Int, Float64}(i => JuMP.dual(c) for (i, c) in conservation)
    mu_e = Dict{Int, Float64}(a => JuMP.dual(c) for (a, c) in loss)
    h = Dict{Int, Float64}(i => JuMP.value(h[i]) for i in h_ids)
    q = Dict{Int, Float64}(a => JuMP.value(q[a]) for a in q_ids)
    qr = Dict{Int, Float64}(i => JuMP.value(qr[i]) for i in ids(wm, n, :reservoirs))
    Delta_v = Dict{Int, Float64}(i => JuMP.value(Delta_v[i]) for i in h_ids)
    Delta_e = Dict{Int, Float64}(a => JuMP.value(Delta_e[a]) for a in q_ids)

    #for i in h_ids
    #    arcs_fr = ref(wm, n, :node_arcs_fr, i)
    #    arcs_to = ref(wm, n, :node_arcs_to, i)
    #    res = ref(wm, n, :node_reservoirs, i)
    #    junc = ref(wm, n, :node_junctions, i)
    #    demands = Dict(k => ref(wm, n, :junctions, k, "demand") for k in junc)

    #    #println(i, " ", arcs_fr, " ", arcs_to)
    #    #println(i, " ", arcs_fr, " ", arcs_to)
    #    ##println(i, " ", arcs_fr, " ", arcs_to)
    #    #println(i, " ", arcs_fr, " ", arcs_to)

    #    lhs = length(arcs_fr) > 0 ? -sum(q[l] for (l, f, t) in arcs_fr) : 0.0
    #    lhs += length(arcs_to) > 0 ? sum(q[l] for (l, f, t) in arcs_to) : 0.0
    #    rhs = length(res) > 0 ? -sum(qr[rid] for rid in res) : 0.0
    #    rhs += length(demands) > 0 ? sum(demand for (jid, demand) in demands) : 0.0

    #    #println(i, " ", lhs, " ", rhs)

    #    #    sum(-qr[rid] for rid in res) +
    #    #    sum(demand for (jid, demand) in demands))
    #end

    # TODO: Find a way to calculate this.
    zeta = compute_zeta(wm, q, mu_e, h, mu_v, lambda_e_p, lambda_e_n, n)
    println(zeta)

    lhs = JuMP.AffExpr(0.0)
    rhs = JuMP.AffExpr(0.0)

    rhs -= zeta * sum(h[i] * junc["demand"] for (i, junc) in ref(wm, n, :junctions))
    rhs += zeta * sum(h[i] * qr[i] for (i, res) in ref(wm, n, :reservoirs))
    rhs -= (1 - zeta) * sum(mu_v[i] * junc["demand"] for (i, junc) in ref(wm, n, :junctions))
    rhs += (1 - zeta) * sum(mu_v[i] * qr[i] for (i, res) in ref(wm, n, :reservoirs))
    rhs += (1 - zeta) * sum(lambda_v_p[i] * JuMP.upper_bound(WMs.var(wm, n, :h, i)) for i in h_ids)
    rhs -= (1 - zeta) * sum(lambda_v_n[i] * JuMP.lower_bound(WMs.var(wm, n, :h, i)) for i in h_ids)

    for a in q_ids
        for (r_id, r) in enumerate(ref(wm, n, :resistance, a))
            i = WMs.ref(wm, n, :links, a)["node_fr"]
            j = WMs.ref(wm, n, :links, a)["node_to"]
            res = WMs.ref(wm, n, :resistance, a)[r_id]

            h_i_ub = JuMP.upper_bound(WMs.var(wm, n, :h, i))
            h_i_lb = JuMP.lower_bound(WMs.var(wm, n, :h, i))
            h_j_ub = JuMP.upper_bound(WMs.var(wm, n, :h, j))
            h_j_lb = JuMP.lower_bound(WMs.var(wm, n, :h, j))

            q_ub = JuMP.upper_bound(WMs.var(wm, n, :qp_ne, a)[r_id])
            q_lb = -JuMP.upper_bound(WMs.var(wm, n, :qn_ne, a)[r_id])

            tau = compute_tau(wm, q, mu_e, h, mu_v, lambda_e_p, lambda_e_n, zeta, a, optimizer, res, n)
            lhs += tau * var(wm, n, :x_res, a)[r_id]
            rhs += (1 - zeta) * var(wm, n, :x_res, a)[r_id] * (lambda_e_p[a] * q_ub)
            rhs -= (1 - zeta) * var(wm, n, :x_res, a)[r_id] * (lambda_e_n[a] * q_lb)

            #rhs += zeta * max(q[a] * (h_i_ub - h_j_lb), q[a] * (h_i_lb - h_j_ub))
            #rhs += zeta * max(q[a] * (h_i_ub - h_j_lb), q[a] * (h_i_lb - h_j_ub))
            #rhs += (1 - zeta) * max(mu_e[a] * (h_i_ub - h_j_lb), mu_e[a] * (h_i_lb - h_j_ub))
        end
    end
   
    c = JuMP.@constraint(wm.model, lhs <= rhs)
    ##println(c)
    ##println(c)

    #xr_ones = Array{JuMP.VariableRef, 1}()
    #xr_zeros = Array{JuMP.VariableRef, 1}()

    #for (a, link) in wm.ref[:nw][n][:links_ne]
    #    xr = JuMP.value.(wm.var[:nw][n][:x_res][a])
    #    r_id = resistance_indices[a]
    #    zero_indices = setdiff(1:length(xr), [r_id])
    #    xr_ones = vcat(xr_ones, wm.var[:nw][n][:x_res][a][r_id])
    #    xr_zeros = vcat(xr_zeros, wm.var[:nw][n][:x_res][a][zero_indices])
    #end

    #no_good_rhs = length(wm.ref[:nw][n][:links_ne]) - 1
    #JuMP.@constraint(wm.model, sum(xr_ones) - sum(xr_zeros) <= no_good_rhs)
end

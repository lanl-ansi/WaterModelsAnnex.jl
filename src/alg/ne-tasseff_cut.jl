function add_tasseff_cut!(wm::GenericWaterModel, optimizer::JuMP.OptimizerFactory, q_sol::Dict{Int, Float64}, h_sol::Dict{Int, Float64}, n::Int=wm.cnw)
    alpha = wm.ref[:nw][n][:alpha]
    resistance_indices = Dict{Int, Int}(a => 0 for a in ids(wm, n, :links_ne))

    for (a, link) in wm.ref[:nw][n][:links_ne]
        x_res, r_id = findmax(JuMP.value.(wm.var[:nw][n][:x_res][a]))
        resistance_indices[a] = r_id
    end

    model = JuMP.Model(JuMP.with_optimizer(optimizer))
    link_ids = collect(WMs.ids(wm, n, :links))
    junction_ids = collect(WMs.ids(wm, n, :junctions))
    reservoir_ids = collect(WMs.ids(wm, n, :reservoirs))
    node_ids = vcat(junction_ids, reservoir_ids)

    fp = JuMP.AffExpr(0.0)
    fd = JuMP.AffExpr(0.0)

    qp_hat = Dict{Int, Array{Float64}}()
    qn_hat = Dict{Int, Array{Float64}}()

    mccormick_1 = Dict{Int, JuMP.ConstraintRef}()
    mccormick_2 = Dict{Int, JuMP.ConstraintRef}()
    mccormick_3 = Dict{Int, JuMP.ConstraintRef}()
    mccormick_4 = Dict{Int, JuMP.ConstraintRef}()

    qp_upper_bound = Dict{Int, JuMP.ConstraintRef}()
    qn_upper_bound = Dict{Int, JuMP.ConstraintRef}()

    gp_definition = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    gn_definition = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()

    gp_definition_rhs = Dict{Int, Array{Float64}}()
    gn_definition_rhs = Dict{Int, Array{Float64}}()

    num_outer_cuts_p = Dict{Int, Int}()
    num_outer_cuts_n = Dict{Int, Int}()

    lambda_p = Dict{Int, Array{JuMP.VariableRef}}()
    lambda_n = Dict{Int, Array{JuMP.VariableRef}}()

    diff_p_definition = Dict{Int, JuMP.ConstraintRef}()
    diff_n_definition = Dict{Int, JuMP.ConstraintRef}()

    lambda_p_definition = Dict{Int, JuMP.ConstraintRef}()
    lambda_n_definition = Dict{Int, JuMP.ConstraintRef}()

    dh = JuMP.@variable(model, [a in link_ids], base_name="dh[$a]")
    qp = JuMP.@variable(model, [a in link_ids], lower_bound=0.0, base_name="qp[$a]")
    qn = JuMP.@variable(model, [a in link_ids], lower_bound=0.0, base_name="qn[$a]")
    gp = JuMP.@variable(model, [a in link_ids], lower_bound=0.0, base_name="gp[$a]")
    gn = JuMP.@variable(model, [a in link_ids], lower_bound=0.0, base_name="gn[$a]")
    qr = JuMP.@variable(model, [i in reservoir_ids], base_name = "qr")
    h = JuMP.@variable(model, [i in node_ids], base_name = "h")

    for (a, link) in wm.ref[:nw][n][:links]
        L = link["length"]
        i = link["node_fr"]
        j = link["node_to"]

        r = resistance_indices[a]
        r_a = wm.ref[:nw][n][:resistance][a][r]

        # Add upper bound constraints.
    	qp_a_r_ub = JuMP.upper_bound(wm.var[:nw][n][:qp_ne][a][r])
    	qp_upper_bound[a] = JuMP.@constraint(model, qp[a] <= qp_a_r_ub)
    	qn_a_r_ub = JuMP.upper_bound(wm.var[:nw][n][:qn_ne][a][r])
        qn_upper_bound[a] = JuMP.@constraint(model, qn[a] <= qn_a_r_ub)

        num_outer_cuts_n[a] = 250
        num_outer_cuts_p[a] = 250

        # Suggest outer-approximation points.
        if q_sol[a] >= 0.0
            first_range = range(0.0, stop = max(0.0, q_sol[a] - 1.0e-3), length = convert(Int, (num_outer_cuts_p[a] - 1) / 2))
            second_range = range(q_sol[a] + 1.0e-3, stop = qp_a_r_ub, length = convert(Int, (num_outer_cuts_p[a] - 1) / 2))
            qp_hat[a] = [first_range; q_sol[a]; second_range]
            qn_hat[a] = range(0.0, stop = qn_a_r_ub, length = num_outer_cuts_n[a])
        else
            first_range = range(0.0, stop = max(0.0, -q_sol[a] - 1.0e-3), length = convert(Int, (num_outer_cuts_n[a] - 1) / 2))
            second_range = range(-q_sol[a] + 1.0e-3, stop = qn_a_r_ub, length = convert(Int, (num_outer_cuts_n[a] - 1) / 2))
            qp_hat[a] = range(0.0, stop = qp_a_r_ub, length = num_outer_cuts_p[a])
            qn_hat[a] = [first_range; -q_sol[a]; second_range]
        end

        #if qp_a_r_ub > 0.0
        #    num_outer_cuts_p[a] = 250
        #    qp_hat[a] = range(0.0, stop=qp_a_r_ub, length=num_outer_cuts_p[a]+2)
        #    qp_hat[a] = qp_hat[a][2:num_outer_cuts_p[a]+1]
        #else
        #    num_outer_cuts_p[a] = 1
        #    qp_hat[a] = [0.0]
        #end

        #if qn_a_r_ub > 0.0
        #    num_outer_cuts_n[a] = 250
        #    qn_hat[a] = range(0.0, stop=qn_a_r_ub, length=num_outer_cuts_n[a]+2)
        #    qn_hat[a] = qn_hat[a][2:num_outer_cuts_n[a]+1]
        #else
        #    num_outer_cuts_n[a] = 1
        #    qn_hat[a] = [0.0]
        #end

        gp_definition[a] = Dict{Int, JuMP.ConstraintRef}()
        gn_definition[a] = Dict{Int, JuMP.ConstraintRef}()
        gp_definition_rhs[a] = zeros(num_outer_cuts_p[a])
        gn_definition_rhs[a] = zeros(num_outer_cuts_n[a])

        # Add positive flow outer-approximation constraints.
        for k in 1:num_outer_cuts_p[a]
            rhs_p = (1.0 - inv(1.0 + alpha)) * qp_hat[a][k]^(1.0 + alpha)
            lhs_p = qp_hat[a][k]^alpha * qp[a] - gp[a]
            gp_definition[a][k] = JuMP.@constraint(model, lhs_p <= rhs_p)
            gp_definition_rhs[a][k] = rhs_p
        end

        # Add negative flow outer-approximation constraints.
        for k in 1:num_outer_cuts_n[a]
            rhs_n = (1.0 - inv(1.0 + alpha)) * qn_hat[a][k]^(1.0 + alpha)
            lhs_n = qn_hat[a][k]^alpha * qn[a] - gn[a]
            gn_definition[a][k] = JuMP.@constraint(model, lhs_n <= rhs_n)
            gn_definition_rhs[a][k] = rhs_n
        end

        fp -= L * r_a * (gp[a] + gn[a])

        lambda_p[a] = JuMP.@variable(model, [k in 1:num_outer_cuts_p[a]], lower_bound=0.0, base_name="lambda_p[$(a)]")
        lambda_n[a] = JuMP.@variable(model, [k in 1:num_outer_cuts_n[a]], lower_bound=0.0, base_name="lambda_n[$(a)]")

        # Add dual constraints of flow variables.
        diff_p_expr = r_a * sum(qp_hat[a].^alpha .* lambda_p[a])
        diff_n_expr = r_a * sum(qn_hat[a].^alpha .* lambda_n[a])

        dh_ub = JuMP.upper_bound(wm.var[:nw][n][:h][i]) - JuMP.lower_bound(wm.var[:nw][n][:h][j])
        dh_lb = JuMP.lower_bound(wm.var[:nw][n][:h][i]) - JuMP.upper_bound(wm.var[:nw][n][:h][j])

        mccormick_1[a] = JuMP.@constraint(model, dh[a] >= dh_lb)
        mccormick_2[a] = JuMP.@constraint(model, dh[a] - (h[i] - h[j]) >= 0.0)
        mccormick_3[a] = JuMP.@constraint(model, dh[a] <= dh_ub)
        mccormick_4[a] = JuMP.@constraint(model, dh[a] - (h[i] - h[j]) <= 0.0)

        diff_p_definition[a] = JuMP.@constraint(model, diff_p_expr - inv(L) * dh[a] >= 0.0)
        diff_n_definition[a] = JuMP.@constraint(model, diff_n_expr + inv(L) * dh[a] >= 0.0)

        # Add dual constraints of g variables.
        lambda_p_definition[a] = JuMP.@constraint(model, sum(lambda_p[a]) <= 1.0)
        lambda_n_definition[a] = JuMP.@constraint(model, sum(lambda_n[a]) <= 1.0)

        fd += L * r_a * sum(lambda_p[a] .* gp_definition_rhs[a])
        fd += L * r_a * sum(lambda_n[a] .* gn_definition_rhs[a])
    end

    # Define flow conservation constraints.
    fp += sum(qr[i] * res["head"] for (i, res) in ref(wm, n, :reservoirs))
    flow_conservation = Dict{Int, JuMP.ConstraintRef}()

    for (i, node) in wm.ref[:nw][n][:nodes]
        junc = ref(wm, n, :node_junctions, i)
        arcs_fr = ref(wm, n, :node_arcs_fr, i)
        arcs_to = ref(wm, n, :node_arcs_to, i)
        res = ref(wm, n, :node_reservoirs, i)
        tank = ref(wm, n, :node_tanks, i)
        demands = Dict(k => ref(wm, n, :junctions, k, "demand") for k in junc)

        flow_conservation[i] = JuMP.@constraint(model,
            sum(-qp[l] + qn[l] for (l, f, t) in arcs_fr) +
            sum( qp[l] - qn[l] for (l, f, t) in arcs_to) ==
            sum(-qr[rid] for rid in res) +
            sum(demand for (jid, demand) in demands))
    end

    h_lower_bound = Dict{Int, JuMP.ConstraintRef}()
    h_upper_bound = Dict{Int, JuMP.ConstraintRef}()

    for (i, node) in wm.ref[:nw][n][:nodes]
        h_i_lb = JuMP.lower_bound(wm.var[:nw][n][:h][i])
        h_i_ub = JuMP.upper_bound(wm.var[:nw][n][:h][i])
        h_lower_bound[i] = JuMP.@constraint(model, h[i] >= h_i_lb)
        h_upper_bound[i] = JuMP.@constraint(model, h[i] <= h_i_ub)
    end

    for (i, junction) in wm.ref[:nw][n][:junctions]
        fd += h[i] * junction["demand"]
    end

    JuMP.@objective(model, MOI.MAX_SENSE, fp - fd)
    strong_duality = JuMP.@constraint(model, fp - fd <= 0.0)
    JuMP.optimize!(model)
    println(JuMP.dual_status(model))

    if JuMP.dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE
        gt_cut_lhs = JuMP.AffExpr(0.0)

        for (a, link) in wm.ref[:nw][n][:links_ne]
            L = link["length"]
            i = link["node_fr"]
            j = link["node_to"]

            r = resistance_indices[a]
            r_a = wm.ref[:nw][n][:resistance][a][r]
            x_r_a = wm.var[:nw][n][:x_res][a][r]

            q_p_a_ub = JuMP.upper_bound(wm.var[:nw][n][:qp_ne][a][r])
            q_n_a_ub = JuMP.upper_bound(wm.var[:nw][n][:qn_ne][a][r])

            dh_ub = JuMP.upper_bound(wm.var[:nw][n][:h][i]) - JuMP.lower_bound(wm.var[:nw][n][:h][j])
            dh_lb = JuMP.lower_bound(wm.var[:nw][n][:h][i]) - JuMP.upper_bound(wm.var[:nw][n][:h][j])

            gt_cut_lhs += JuMP.dual(mccormick_1[a]) * (dh_lb * x_r_a)
            gt_cut_lhs += JuMP.dual(mccormick_2[a]) * (dh_ub * x_r_a - dh_ub)
            gt_cut_lhs += JuMP.dual(mccormick_3[a]) * (dh_ub * x_r_a)
            gt_cut_lhs += JuMP.dual(mccormick_4[a]) * (dh_lb * x_r_a - dh_lb)
            gt_cut_lhs += JuMP.dual(lambda_p_definition[a]) * x_r_a
            gt_cut_lhs += JuMP.dual(lambda_n_definition[a]) * x_r_a

            for k in 1:num_outer_cuts_p[a]
                gt_cut_lhs += JuMP.dual(gp_definition[a][k]) * gp_definition_rhs[a][k]
            end

            for k in 1:num_outer_cuts_n[a]
                gt_cut_lhs += JuMP.dual(gn_definition[a][k]) * gn_definition_rhs[a][k]
            end

            gt_cut_lhs += JuMP.dual(qp_upper_bound[a]) * q_p_a_ub * x_r_a
            gt_cut_lhs += JuMP.dual(qn_upper_bound[a]) * q_n_a_ub * x_r_a
        end

        for (i, nodes) in wm.ref[:nw][n][:nodes]
            gt_cut_lhs += JuMP.dual(h_lower_bound[i]) * JuMP.lower_bound(wm.var[:nw][n][:h][i])
            gt_cut_lhs += JuMP.dual(h_upper_bound[i]) * JuMP.upper_bound(wm.var[:nw][n][:h][i])

            if i in keys(wm.ref[:nw][n][:junctions])
                junction = wm.ref[:nw][n][:junctions][i]
                gt_cut_lhs += JuMP.dual(flow_conservation[i]) * junction["demand"]
            end
        end

        #add_no_good_cut!(wm, optimizer, n)
        c = JuMP.@constraint(wm.model, gt_cut_lhs >= 0.0)
        #println(c)
    else
        #add_no_good_cut!(wm, optimizer, n)
    end
end

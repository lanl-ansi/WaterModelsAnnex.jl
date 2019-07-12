function add_tasseff_cut!(wm::GenericWaterModel, optimizer::JuMP.OptimizerFactory, n::Int=wm.cnw)
    alpha = wm.ref[:nw][n][:alpha]
    resistance_indices = Dict{Int, Int}(a => 0 for a in ids(wm, n, :links_ne))

    for (a, link) in wm.ref[:nw][n][:links_ne]
        x_res, r_id = findmax(JuMP.value.(wm.var[:nw][wm.cnw][:x_res][a]))
        resistance_indices[a] = r_id
    end

    model = JuMP.Model(JuMP.with_optimizer(optimizer))

    fp = JuMP.AffExpr(0.0)
    fd = JuMP.AffExpr(0.0)

    qp = Dict{Int, Array{JuMP.VariableRef}}()
    qn = Dict{Int, Array{JuMP.VariableRef}}()

    gp = Dict{Int, Array{JuMP.VariableRef}}()
    gn = Dict{Int, Array{JuMP.VariableRef}}()

    qp_hat = Dict{Int, Dict{Int, Array{Float64}}}()
    qn_hat = Dict{Int, Dict{Int, Array{Float64}}}()

    mccormick_1 = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    mccormick_2 = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    mccormick_3 = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    mccormick_4 = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()

    qp_upper_bound = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    qn_upper_bound = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()

    gp_definition = Dict{Int, Dict{Int, Dict{Int, JuMP.ConstraintRef}}}()
    gn_definition = Dict{Int, Dict{Int, Dict{Int, JuMP.ConstraintRef}}}()

    gp_definition_rhs = Dict{Int, Dict{Int, Array{Float64}}}()
    gn_definition_rhs = Dict{Int, Dict{Int, Array{Float64}}}()

    num_outer_cuts_p = Dict{Int, Dict{Int, Int}}()
    num_outer_cuts_n = Dict{Int, Dict{Int, Int}}()

    lambda_p = Dict{Int, Dict{Int, Array{JuMP.VariableRef}}}()
    lambda_n = Dict{Int, Dict{Int, Array{JuMP.VariableRef}}}()

    diff_p_definition = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    diff_n_definition = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()

    lambda_p_definition = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    lambda_n_definition = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()

    junction_ids = WMs.ids(wm, n, :junctions)
    reservoir_ids = WMs.ids(wm, n, :reservoirs)
    nodes = [junction_ids; reservoir_ids]
    h = JuMP.@variable(model, [i in junction_ids], base_name = "h")
    dh = Dict{Int, Array{JuMP.VariableRef}}()

    for (a, link) in wm.ref[:nw][n][:links]
        L = link["length"]
        i = link["f_id"]
        j = link["t_id"]

        num_resistances = length(wm.ref[:nw][n][:resistance][a])

        dh[a] = JuMP.@variable(model, [r in 1:num_resistances], base_name = "dh[$a]")

        qp[a] = JuMP.@variable(model, [r in 1:num_resistances], lower_bound = 0.0, base_name = "qp[$a]")
        qn[a] = JuMP.@variable(model, [r in 1:num_resistances], lower_bound = 0.0, base_name = "qn[$a]")

        gp[a] = JuMP.@variable(model, [r in 1:num_resistances], lower_bound = 0.0, base_name = "gp[$a]")
        gn[a] = JuMP.@variable(model, [r in 1:num_resistances], lower_bound = 0.0, base_name = "gn[$a]")

        gp_definition[a] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        gn_definition[a] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()

        gp_definition_rhs[a] = Dict{Int, Array{Float64}}()
        gn_definition_rhs[a] = Dict{Int, Array{Float64}}()

        qp_upper_bound[a] = Dict{Int, JuMP.ConstraintRef}()
        qn_upper_bound[a] = Dict{Int, JuMP.ConstraintRef}()

        mccormick_1[a] = Dict{Int, JuMP.ConstraintRef}()
        mccormick_2[a] = Dict{Int, JuMP.ConstraintRef}()
        mccormick_3[a] = Dict{Int, JuMP.ConstraintRef}()
        mccormick_4[a] = Dict{Int, JuMP.ConstraintRef}()

        qp_hat[a] = Dict{Int, Array{Float64}}()
        qn_hat[a] = Dict{Int, Array{Float64}}()

        num_outer_cuts_p[a] = Dict{Int, Int}()
        num_outer_cuts_n[a] = Dict{Int, Int}()

        lambda_p[a] = Dict{Int, Array{JuMP.VariableRef}}()
        lambda_n[a] = Dict{Int, Array{JuMP.VariableRef}}()

        lambda_p_definition[a] = Dict{Int, JuMP.ConstraintRef}()
        lambda_n_definition[a] = Dict{Int, JuMP.ConstraintRef}()

        diff_p_definition[a] = Dict{Int, JuMP.ConstraintRef}()
        diff_n_definition[a] = Dict{Int, JuMP.ConstraintRef}()

        for r in 1:length(wm.ref[:nw][n][:resistance][a])
            x_sel = r == resistance_indices[a] ? 1.0 : 0.0
            r_a = wm.ref[:nw][n][:resistance][a][r]

            # Add upper bound constraints.
    		qp_a_r_ub = JuMP.upper_bound(wm.var[:nw][n][:qp_ne][a][r])
    		qp_upper_bound[a][r] = JuMP.@constraint(model, qp[a][r] <= qp_a_r_ub * x_sel)
    		qn_a_r_ub = JuMP.upper_bound(wm.var[:nw][n][:qn_ne][a][r])
            qn_upper_bound[a][r] = JuMP.@constraint(model, qn[a][r] <= qn_a_r_ub * x_sel)

            # Suggest outer-approximation points.
            if qp_a_r_ub > 0.0
                num_outer_cuts_p[a][r] = 25
                qp_hat[a][r] = range(0.0, stop=qp_a_r_ub, length=num_outer_cuts_p[a][r]+2)
                qp_hat[a][r] = qp_hat[a][r][2:num_outer_cuts_p[a][r]+1]
            else
                num_outer_cuts_p[a][r] = 1
                qp_hat[a][r] = [0.0]
            end

            if qn_a_r_ub > 0.0
                num_outer_cuts_n[a][r] = 25
                qn_hat[a][r] = range(0.0, stop=qn_a_r_ub, length=num_outer_cuts_n[a][r]+2)
                qn_hat[a][r] = qn_hat[a][r][2:num_outer_cuts_n[a][r]+1]
            else
                num_outer_cuts_n[a][r] = 1
                qn_hat[a][r] = [0.0]
            end

    		gp_definition[a][r] = Dict{Int, JuMP.ConstraintRef}()
    		gn_definition[a][r] = Dict{Int, JuMP.ConstraintRef}()

            gp_definition_rhs[a][r] = zeros(num_outer_cuts_p[a][r])
            gn_definition_rhs[a][r] = zeros(num_outer_cuts_n[a][r])

            # Add positive flow outer-approximation constraints.
            for k in 1:num_outer_cuts_p[a][r]
                rhs_p = (1.0 - inv(1.0 + alpha)) * qp_hat[a][r][k]^(1.0 + alpha)
                lhs_p = qp_hat[a][r][k]^alpha * qp[a][r] - gp[a][r]
                gp_definition[a][r][k] = JuMP.@constraint(model, lhs_p <= rhs_p)
                gp_definition_rhs[a][r][k] = rhs_p
            end

            # Add negative flow outer-approximation constraints.
            for k in 1:num_outer_cuts_n[a][r]
                rhs_n = (1.0 - inv(1.0 + alpha)) * qn_hat[a][r][k]^(1.0 + alpha)
                lhs_n = qn_hat[a][r][k]^alpha * qn[a][r] - gn[a][r]
                gn_definition[a][r][k] = JuMP.@constraint(model, lhs_n <= rhs_n)
                gn_definition_rhs[a][r][k] = rhs_n
            end

            fp -= L * wm.ref[:nw][n][:resistance][a][r] * (gp[a][r] + gn[a][r])

            lambda_p[a][r] = JuMP.@variable(model, [k in 1:num_outer_cuts_p[a][r]], lower_bound = 0.0, base_name = "lambda_p[$(a)][$(r)]")
            lambda_n[a][r] = JuMP.@variable(model, [k in 1:num_outer_cuts_n[a][r]], lower_bound = 0.0, base_name = "lambda_n[$(a)][$(r)]")

            # Add dual constraints of flow variables.
            diff_p_expr = r_a * sum(qp_hat[a][r].^alpha .* lambda_p[a][r])
            diff_n_expr = r_a * sum(qn_hat[a][r].^alpha .* lambda_n[a][r])

            if i in reservoir_ids
                h_i_res = wm.ref[:nw][n][:reservoirs][i]["head"]

                dh_ub = h_i_res - JuMP.lower_bound(wm.var[:nw][n][:h][j])
                dh_lb = h_i_res - JuMP.upper_bound(wm.var[:nw][n][:h][j])

                mccormick_1[a][r] = JuMP.@constraint(model, dh[a][r] >= dh_lb * x_sel)
                mccormick_2[a][r] = JuMP.@constraint(model, dh[a][r] + h[j] >= dh_ub * x_sel + h_i_res - dh_ub)
                mccormick_3[a][r] = JuMP.@constraint(model, dh[a][r] <= x_sel * dh_ub)
                mccormick_4[a][r] = JuMP.@constraint(model, dh[a][r] + h[j] <= h_i_res + dh_lb * x_sel - dh_lb)

                diff_p_definition[a][r] = JuMP.@constraint(model, diff_p_expr - inv(L) * dh[a][r] >= 0.0)
                diff_n_definition[a][r] = JuMP.@constraint(model, diff_n_expr + inv(L) * dh[a][r] >= 0.0)
            elseif i in junction_ids
                dh_ub = JuMP.upper_bound(wm.var[:nw][n][:h][i]) - JuMP.lower_bound(wm.var[:nw][n][:h][j])
                dh_lb = JuMP.lower_bound(wm.var[:nw][n][:h][i]) - JuMP.upper_bound(wm.var[:nw][n][:h][j])

                mccormick_1[a][r] = JuMP.@constraint(model, dh[a][r] >= dh_lb * x_sel)
                mccormick_2[a][r] = JuMP.@constraint(model, dh[a][r] - (h[i] - h[j]) >= dh_ub * x_sel - dh_ub)
                mccormick_3[a][r] = JuMP.@constraint(model, dh[a][r] <= dh_ub * x_sel)
                mccormick_4[a][r] = JuMP.@constraint(model, dh[a][r] - (h[i] - h[j]) <= dh_lb * x_sel - dh_lb)

                diff_p_definition[a][r] = JuMP.@constraint(model, diff_p_expr - inv(L) * dh[a][r] >= 0.0)
                diff_n_definition[a][r] = JuMP.@constraint(model, diff_n_expr + inv(L) * dh[a][r] >= 0.0)
            end

            # Add dual constraints of g variables.
            lambda_p_definition[a][r] = JuMP.@constraint(model, sum(lambda_p[a][r]) <= x_sel)
            lambda_n_definition[a][r] = JuMP.@constraint(model, sum(lambda_n[a][r]) <= x_sel)

            fd += L * r_a * sum(lambda_p[a][r] .* gp_definition_rhs[a][r])
            fd += L * r_a * sum(lambda_n[a][r] .* gn_definition_rhs[a][r])
        end
    end

    for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
        for (a, link) in filter(a -> i == a.second["f_id"], wm.ref[:nw][n][:links])
            r_id = resistance_indices[a]
            fp += reservoir["head"] * (qp[a][r_id] - qn[a][r_id])
        end

        for (a, link) in filter(a -> i == a.second["t_id"], wm.ref[:nw][n][:links])
            r_id = resistance_indices[a]
            fp -= reservoir["head"] * (qp[a][r_id] - qn[a][r_id])
        end
    end

    # Define flow conservation constraints.
    flow_conservation = Dict{Int, JuMP.ConstraintRef}()

    for (i, junction) in wm.ref[:nw][n][:junctions]
        flow_sum = JuMP.AffExpr(0.0)

        for (a, conn) in filter(a -> i == a.second["t_id"], wm.ref[:nw][n][:links])
            r_id = resistance_indices[a]
            flow_sum += qp[a][r_id] - qn[a][r_id]
        end

        for (a, conn) in filter(a -> i == a.second["f_id"], wm.ref[:nw][n][:links])
            r_id = resistance_indices[a]
            flow_sum -= qp[a][r_id] - qn[a][r_id]
        end

        flow_conservation[i] = JuMP.@constraint(model, flow_sum == junction["demand"])
    end

    h_lower_bound = Dict{Int, JuMP.ConstraintRef}()
    h_upper_bound = Dict{Int, JuMP.ConstraintRef}()

    for (i, junction) in wm.ref[:nw][n][:junctions]
        fd += h[i] * junction["demand"]
        h_i_lb = JuMP.lower_bound(wm.var[:nw][n][:h][i])
        h_i_ub = JuMP.upper_bound(wm.var[:nw][n][:h][i])
        h_lower_bound[i] = JuMP.@constraint(model, h[i] >= h_i_lb)
        h_upper_bound[i] = JuMP.@constraint(model, h[i] <= h_i_ub)
    end

    JuMP.@objective(model, MOI.MAX_SENSE, fp - fd)
    strong_duality = JuMP.@constraint(model, fp - fd == 0.0)
    JuMP.optimize!(model)

    if JuMP.dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE
        gt_cut_lhs = JuMP.AffExpr(0.0)

        for (a, link) in wm.ref[:nw][n][:links_ne]
            L = link["length"]
            i = link["f_id"]
            j = link["t_id"]

            for r in 1:length(wm.ref[:nw][n][:resistance][a])
                r_a = wm.ref[:nw][n][:resistance][a][r]
                x_r_a = wm.var[:nw][n][:x_res][a][r]
                x_sel = r == resistance_indices[a] ? 1.0 : 0.0

                q_p_a_ub = JuMP.upper_bound(wm.var[:nw][n][:qp_ne][a][r])
                q_n_a_ub = JuMP.upper_bound(wm.var[:nw][n][:qn_ne][a][r])

                if i in collect(reservoir_ids)
                    reservoir = wm.ref[:nw][n][:reservoirs][i]
                    h_i_res = reservoir["head"]

                    dh_ub = h_i_res - JuMP.lower_bound(wm.var[:nw][n][:h][j])
                    dh_lb = h_i_res - JuMP.upper_bound(wm.var[:nw][n][:h][j])

                    gt_cut_lhs += JuMP.dual(mccormick_1[a][r]) * (dh_lb * x_r_a)
                    gt_cut_lhs += JuMP.dual(mccormick_2[a][r]) * (dh_ub * x_r_a + h_i_res - dh_ub)
                    gt_cut_lhs += JuMP.dual(mccormick_3[a][r]) * (dh_ub * x_r_a)
                    gt_cut_lhs += JuMP.dual(mccormick_4[a][r]) * (dh_lb * x_r_a + h_i_res - dh_lb)

                    #gt_cut_lhs -= JuMP.dual(diff_p_definition[a][r]) * inv(L) * x_sel * h_i_res
                    #gt_cut_lhs += JuMP.dual(diff_n_definition[a][r]) * inv(L) * x_sel * h_i_res
                else
                    dh_ub = JuMP.upper_bound(wm.var[:nw][n][:h][i]) - JuMP.lower_bound(wm.var[:nw][n][:h][j])
                    dh_lb = JuMP.lower_bound(wm.var[:nw][n][:h][i]) - JuMP.upper_bound(wm.var[:nw][n][:h][j])

                    gt_cut_lhs += JuMP.dual(mccormick_1[a][r]) * (dh_lb * x_r_a)
                    gt_cut_lhs += JuMP.dual(mccormick_2[a][r]) * (dh_ub * x_r_a - dh_ub)
                    gt_cut_lhs += JuMP.dual(mccormick_3[a][r]) * (dh_ub * x_r_a)
                    gt_cut_lhs += JuMP.dual(mccormick_4[a][r]) * (dh_lb * x_r_a - dh_lb)
                end

                gt_cut_lhs += JuMP.dual(lambda_p_definition[a][r]) * x_r_a
                gt_cut_lhs += JuMP.dual(lambda_n_definition[a][r]) * x_r_a

                for k in 1:num_outer_cuts_p[a][r]
                    gt_cut_lhs += JuMP.dual(gp_definition[a][r][k]) * gp_definition_rhs[a][r][k]
                end

                for k in 1:num_outer_cuts_n[a][r]
                    gt_cut_lhs += JuMP.dual(gn_definition[a][r][k]) * gn_definition_rhs[a][r][k]
                end

                gt_cut_lhs += JuMP.dual(qp_upper_bound[a][r]) * q_p_a_ub * x_r_a
                gt_cut_lhs += JuMP.dual(qn_upper_bound[a][r]) * q_n_a_ub * x_r_a
            end
        end

        for (i, junction) in wm.ref[:nw][n][:junctions]
            gt_cut_lhs += JuMP.dual(h_lower_bound[i]) * JuMP.lower_bound(wm.var[:nw][n][:h][i])
            gt_cut_lhs += JuMP.dual(h_upper_bound[i]) * JuMP.upper_bound(wm.var[:nw][n][:h][i])
            gt_cut_lhs += JuMP.dual(flow_conservation[i]) * junction["demand"]
        end

        add_no_good_cut!(wm, optimizer, n)
        JuMP.@constraint(wm.model, gt_cut_lhs >= 0.0)
    else
        add_no_good_cut!(wm, optimizer, n)
    end
end

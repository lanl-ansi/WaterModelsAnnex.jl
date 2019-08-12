using Ipopt

function solve_without_bounds(wm::WMs.GenericWaterModel, optimizer, n::Int=wm.cnw)
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
        i = link["node_fr"]
        j = link["node_to"]

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
                qp_hat[a][r] = range(0.0, stop=qp_a_r_ub, length=num_outer_cuts_p[a][r])
                #qp_hat[a][r] = qp_hat[a][r][2:num_outer_cuts_p[a][r]+1]
            else
                num_outer_cuts_p[a][r] = 1
                qp_hat[a][r] = [0.0]
            end

            if qn_a_r_ub > 0.0
                num_outer_cuts_n[a][r] = 25
                qn_hat[a][r] = range(0.0, stop=qn_a_r_ub, length=num_outer_cuts_n[a][r])
                #qn_hat[a][r] = qn_hat[a][r][2:num_outer_cuts_n[a][r]+1]
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
                gp_definition_rhs[a][r][k] = rhs_p * x_sel
            end

            # Add negative flow outer-approximation constraints.
            for k in 1:num_outer_cuts_n[a][r]
                rhs_n = (1.0 - inv(1.0 + alpha)) * qn_hat[a][r][k]^(1.0 + alpha)
                lhs_n = qn_hat[a][r][k]^alpha * qn[a][r] - gn[a][r]
                gn_definition[a][r][k] = JuMP.@constraint(model, lhs_n <= rhs_n)
                gn_definition_rhs[a][r][k] = rhs_n * x_sel
            end

            fp -= L * wm.ref[:nw][n][:resistance][a][r] * (gp[a][r] + gn[a][r])

            lambda_p[a][r] = JuMP.@variable(model, [k in 1:num_outer_cuts_p[a][r]], lower_bound = 0.0, base_name = "lambda_p[$(a)][$(r)]")
            lambda_n[a][r] = JuMP.@variable(model, [k in 1:num_outer_cuts_n[a][r]], lower_bound = 0.0, base_name = "lambda_n[$(a)][$(r)]")

            # Add dual constraints of flow variables.
            diff_p_expr = r_a * sum(qp_hat[a][r].^alpha .* lambda_p[a][r])
            diff_n_expr = r_a * sum(qn_hat[a][r].^alpha .* lambda_n[a][r])

            if i in reservoir_ids
                h_i_res = wm.ref[:nw][n][:reservoirs][i]["head"]

                dh_ub = 5 * (h_i_res - JuMP.lower_bound(wm.var[:nw][n][:h][j]))
                dh_lb = 5 * (h_i_res - JuMP.upper_bound(wm.var[:nw][n][:h][j]))

                mccormick_1[a][r] = JuMP.@constraint(model, dh[a][r] >= dh_lb * x_sel)
                mccormick_2[a][r] = JuMP.@constraint(model, dh[a][r] + h[j] >= dh_ub * x_sel + h_i_res - dh_ub)
                mccormick_3[a][r] = JuMP.@constraint(model, dh[a][r] <= x_sel * dh_ub)
                mccormick_4[a][r] = JuMP.@constraint(model, dh[a][r] + h[j] <= h_i_res + dh_lb * x_sel - dh_lb)

                diff_p_definition[a][r] = JuMP.@constraint(model, diff_p_expr - inv(L) * dh[a][r] >= 0.0)
                diff_n_definition[a][r] = JuMP.@constraint(model, diff_n_expr + inv(L) * dh[a][r] >= 0.0)
            elseif i in junction_ids
                dh_ub = 5 * (JuMP.upper_bound(wm.var[:nw][n][:h][i]) - JuMP.lower_bound(wm.var[:nw][n][:h][j]))
                dh_lb = 5 * (JuMP.lower_bound(wm.var[:nw][n][:h][i]) - JuMP.upper_bound(wm.var[:nw][n][:h][j]))

                mccormick_1[a][r] = JuMP.@constraint(model, dh[a][r] >= dh_lb * x_sel)
                mccormick_2[a][r] = JuMP.@constraint(model, dh[a][r] - (h[i] - h[j]) >= dh_ub * x_sel - dh_ub)
                mccormick_3[a][r] = JuMP.@constraint(model, dh[a][r] <= dh_ub * x_sel)
                mccormick_4[a][r] = JuMP.@constraint(model, dh[a][r] - (h[i] - h[j]) <= dh_lb * x_sel - dh_lb)

                diff_p_definition[a][r] = JuMP.@constraint(model, diff_p_expr - inv(L) * dh[a][r] >= 0.0)
                diff_n_definition[a][r] = JuMP.@constraint(model, diff_n_expr + inv(L) * dh[a][r] >= 0.0)
            end

            # Add dual constraints of g variables.
            lambda_p_definition[a][r] = JuMP.@constraint(model, sum(lambda_p[a][r]) == x_sel)
            lambda_n_definition[a][r] = JuMP.@constraint(model, sum(lambda_n[a][r]) == x_sel)

            fd += L * r_a * sum(lambda_p[a][r] .* gp_definition_rhs[a][r])
            fd += L * r_a * sum(lambda_n[a][r] .* gn_definition_rhs[a][r])
        end
    end

    fp += sum(wm.var[:nw][n][:qr][i] * res["head"] for (i, res) in ref(wm, n, :reservoirs))

    #for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
    #    for (a, link) in filter(a -> i == a.second["node_fr"], wm.ref[:nw][n][:links])
    #        r_id = resistance_indices[a]
    #        fp += reservoir["head"] * (qp[a][r_id] - qn[a][r_id])
    #    end

    #    for (a, link) in filter(a -> i == a.second["node_to"], wm.ref[:nw][n][:links])
    #        r_id = resistance_indices[a]
    #        fp -= reservoir["head"] * (qp[a][r_id] - qn[a][r_id])
    #    end
    #end

    # Define flow conservation constraints.
    flow_conservation = Dict{Int, JuMP.ConstraintRef}()

    for (i, junction) in wm.ref[:nw][n][:junctions]
        flow_sum = JuMP.AffExpr(0.0)

        for (a, conn) in filter(a -> i == a.second["node_to"], wm.ref[:nw][n][:links])
            r_id = resistance_indices[a]
            flow_sum += qp[a][r_id] - qn[a][r_id]
        end

        for (a, conn) in filter(a -> i == a.second["node_fr"], wm.ref[:nw][n][:links])
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
        #h_lower_bound[i] = JuMP.@constraint(model, h[i] >= h_i_lb)
        #h_upper_bound[i] = JuMP.@constraint(model, h[i] <= h_i_ub)
    end

    JuMP.@objective(model, MOI.MIN_SENSE, fd)
    #strong_duality = JuMP.@constraint(model, fd - fp == 0.0)
    JuMP.optimize!(model)

    println("HEAD SOLUTIONS FROM LP")
    for (i, junction) in wm.ref[:nw][n][:junctions]
        h_i_lb = JuMP.lower_bound(wm.var[:nw][n][:h][i])
        h_i_ub = JuMP.upper_bound(wm.var[:nw][n][:h][i])
        println(i, " ", h_i_lb, " ", h_i_ub, " ", JuMP.value(h[i]))
    end

    println("")

    println("OBJECTIVE VALUES FROM LP")
    println(JuMP.value(fp), " ", JuMP.value(fd))

    println("")

    println("LAMBDA SOLUTIONS")
    for (a, link) in wm.ref[:nw][wm.cnw][:links]
        for r in 1:length(wm.ref[:nw][wm.cnw][:resistance][a])
            x_sel = r == resistance_indices[a] ? 1.0 : 0.0

            if x_sel == 1.0
                r_a = wm.ref[:nw][wm.cnw][:resistance][a][r]
                println(a, " + ", JuMP.value.(lambda_p[a][r]), " ", sum(JuMP.value.(lambda_p[a][r])))
                println(a, " - ", JuMP.value.(lambda_n[a][r]), " ", sum(JuMP.value.(lambda_n[a][r])))
            end
        end
    end

    println("")
end

function add_system!(wm::WMs.GenericWaterModel)
    n = wm.cnw
    link_ids = collect(ids(wm, n, :links_ne))
    alpha = wm.ref[:nw][n][:alpha]
    num_links = length(wm.ref[:nw][n][:links_ne])

    h = wm.var[:nw][n][:h]
    qp = wm.var[:nw][n][:qp_ne]
    qn = wm.var[:nw][n][:qn_ne]
    dhp = wm.var[:nw][n][:dhp]
    dhn = wm.var[:nw][n][:dhn]

    wm.var[:nw][n][:fd] = JuMP.AffExpr(0.0)
    wm.var[:nw][n][:fp] = JuMP.AffExpr(0.0)

    gp = Dict{Int, Array{JuMP.VariableRef}}()
    gn = Dict{Int, Array{JuMP.VariableRef}}()

    qp_hat = Dict{Int, Dict{Int, Array{Float64}}}()
    qn_hat = Dict{Int, Dict{Int, Array{Float64}}}()

    mccormick_1 = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    mccormick_2 = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    mccormick_3 = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    mccormick_4 = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()

    gp_definition = Dict{Int, Dict{Int, Dict{Int, JuMP.ConstraintRef}}}()
    gn_definition = Dict{Int, Dict{Int, Dict{Int, JuMP.ConstraintRef}}}()

    gp_definition_rhs = Dict{Int, Dict{Int, Array{Float64}}}()
    gn_definition_rhs = Dict{Int, Dict{Int, Array{Float64}}}()

    num_outer_cuts_p = Dict{Int, Dict{Int, Int}}()
    num_outer_cuts_n = Dict{Int, Dict{Int, Int}}()

    wm.var[:nw][n][:lambda_p] = Dict{Int, Dict{Int, Array{JuMP.VariableRef}}}()
    wm.var[:nw][n][:lambda_n] = Dict{Int, Dict{Int, Array{JuMP.VariableRef}}}()

    diff_p_definition = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    diff_n_definition = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()

    lambda_p_definition = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
    lambda_n_definition = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()

    junction_ids = WMs.ids(wm, n, :junctions)
    reservoir_ids = WMs.ids(wm, n, :reservoirs)
    nodes = [junction_ids; reservoir_ids]
    #dh = Dict{Int, Array{JuMP.VariableRef}}()
    
    dhp_xr = Dict{Int, Array{JuMP.VariableRef}}()
    dhn_xr = Dict{Int, Array{JuMP.VariableRef}}()

    for (a, link) in wm.ref[:nw][n][:links]
        L = link["length"]
        i = link["node_fr"]
        j = link["node_to"]

        dhp_ub = JuMP.upper_bound(dhp[a])
        dhn_ub = JuMP.upper_bound(dhn[a])

        x_dir = wm.var[:nw][n][:x_dir][a]
        num_resistances = length(wm.ref[:nw][n][:resistance][a])

        dhp_xr[a] = JuMP.@variable(wm.model, [r in 1:num_resistances], base_name = "dhp_xr[$a]", start = 0.0, lower_bound = 0.0)
        dhn_xr[a] = JuMP.@variable(wm.model, [r in 1:num_resistances], base_name = "dhn_xr[$a]", start = 0.0, lower_bound = 0.0)

        #dh[a] = JuMP.@variable(wm.model, [r in 1:num_resistances], base_name = "dh[$a]", start = 0.0)
        gp[a] = JuMP.@variable(wm.model, [r in 1:num_resistances], lower_bound = 0.0, base_name = "gp[$a]", start = 0.0)
        gn[a] = JuMP.@variable(wm.model, [r in 1:num_resistances], lower_bound = 0.0, base_name = "gn[$a]", start = 0.0)

        gp_definition[a] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()
        gn_definition[a] = Dict{Int, Dict{Int, JuMP.ConstraintRef}}()

        gp_definition_rhs[a] = Dict{Int, Array{Float64}}()
        gn_definition_rhs[a] = Dict{Int, Array{Float64}}()

        mccormick_1[a] = Dict{Int, JuMP.ConstraintRef}()
        mccormick_2[a] = Dict{Int, JuMP.ConstraintRef}()
        mccormick_3[a] = Dict{Int, JuMP.ConstraintRef}()
        mccormick_4[a] = Dict{Int, JuMP.ConstraintRef}()

        qp_hat[a] = Dict{Int, Array{Float64}}()
        qn_hat[a] = Dict{Int, Array{Float64}}()

        num_outer_cuts_p[a] = Dict{Int, Int}()
        num_outer_cuts_n[a] = Dict{Int, Int}()

        wm.var[:nw][n][:lambda_p][a] = Dict{Int, Array{JuMP.VariableRef}}()
        wm.var[:nw][n][:lambda_n][a] = Dict{Int, Array{JuMP.VariableRef}}()

        lambda_p_definition[a] = Dict{Int, JuMP.ConstraintRef}()
        lambda_n_definition[a] = Dict{Int, JuMP.ConstraintRef}()

        diff_p_definition[a] = Dict{Int, JuMP.ConstraintRef}()
        diff_n_definition[a] = Dict{Int, JuMP.ConstraintRef}()

        for r in 1:length(wm.ref[:nw][n][:resistance][a])
            r_a = wm.ref[:nw][n][:resistance][a][r]
            x_r_a = wm.var[:nw][n][:x_res][a][r]

            # Add upper bound constraints.
    		qp_a_r_ub = JuMP.upper_bound(wm.var[:nw][n][:qp_ne][a][r])
    		qn_a_r_ub = JuMP.upper_bound(wm.var[:nw][n][:qn_ne][a][r])

            # Suggest outer-approximation points.
            if qp_a_r_ub > 0.0
                num_outer_cuts_p[a][r] = 100
                qp_hat[a][r] = range(0.0, stop=qp_a_r_ub, length=num_outer_cuts_p[a][r])
                #qp_hat[a][r] = range(0.0, stop=qp_a_r_ub, length=num_outer_cuts_p[a][r]+2)
                #qp_hat[a][r] = qp_hat[a][r][2:num_outer_cuts_p[a][r]+1]
            else
                num_outer_cuts_p[a][r] = 1
                qp_hat[a][r] = [0.0]
            end

            if qn_a_r_ub > 0.0
                num_outer_cuts_n[a][r] = 100
                qn_hat[a][r] = range(0.0, stop=qn_a_r_ub, length=num_outer_cuts_n[a][r])
                #qn_hat[a][r] = range(0.0, stop=qn_a_r_ub, length=num_outer_cuts_n[a][r]+2)
                #qn_hat[a][r] = qn_hat[a][r][2:num_outer_cuts_n[a][r]+1]
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
                gp_definition[a][r][k] = JuMP.@constraint(wm.model, lhs_p <= rhs_p)
                gp_definition_rhs[a][r][k] = rhs_p
                JuMP.@constraint(wm.model, gp[a][r] <= inv(1.0 + alpha) * qp_a_r_ub^(1.0 + alpha) * x_r_a)
                JuMP.@constraint(wm.model, gp[a][r] <= inv(1.0 + alpha) * qp_a_r_ub^(1.0 + alpha) * x_dir)
            end

            # Add negative flow outer-approximation constraints.
            for k in 1:num_outer_cuts_n[a][r]
                rhs_n = (1.0 - inv(1.0 + alpha)) * qn_hat[a][r][k]^(1.0 + alpha)
                lhs_n = qn_hat[a][r][k]^alpha * qn[a][r] - gn[a][r]
                gn_definition[a][r][k] = JuMP.@constraint(wm.model, lhs_n <= rhs_n)
                gn_definition_rhs[a][r][k] = rhs_n
                JuMP.@constraint(wm.model, gn[a][r] <= inv(1.0 + alpha) * qn_a_r_ub^(1.0 + alpha) * x_r_a)
                JuMP.@constraint(wm.model, gn[a][r] <= inv(1.0 + alpha) * qn_a_r_ub^(1.0 + alpha) * (1 - x_dir))
            end

            wm.var[:nw][n][:fp] -= L * r_a * (gp[a][r] + gn[a][r])

            wm.var[:nw][n][:lambda_p][a][r] = JuMP.@variable(wm.model, [k in 1:num_outer_cuts_p[a][r]], lower_bound=0.0,
                                                             upper_bound = 1.0, base_name="lambda_p[$(a)][$(r)]", start=0.0)
            wm.var[:nw][n][:lambda_n][a][r] = JuMP.@variable(wm.model, [k in 1:num_outer_cuts_n[a][r]], lower_bound=0.0,
                                                             upper_bound = 1.0, base_name="lambda_n[$(a)][$(r)]", start=0.0)

            # Add dual constraints of flow variables.
            diff_p_expr = r_a .* sum(qp_hat[a][r].^alpha .* wm.var[:nw][n][:lambda_p][a][r])
            diff_n_expr = r_a .* sum(qn_hat[a][r].^alpha .* wm.var[:nw][n][:lambda_n][a][r])

            JuMP.@constraint(wm.model, dhp_xr[a][r] <= dhp_ub * x_r_a)
            JuMP.@constraint(wm.model, dhp_xr[a][r] <= dhp_ub * x_dir)
            JuMP.@constraint(wm.model, dhn_xr[a][r] <= dhn_ub * x_r_a)
            JuMP.@constraint(wm.model, dhn_xr[a][r] <= dhn_ub * (1 - x_dir))
            diff_p_definition[a][r] = JuMP.@constraint(wm.model, diff_p_expr - inv(L) * dhp_xr[a][r] >= 0.0)
            diff_n_definition[a][r] = JuMP.@constraint(wm.model, diff_n_expr - inv(L) * dhn_xr[a][r] >= 0.0)

            #if i in reservoir_ids
            #    #h_i_res = wm.ref[:nw][n][:reservoirs][i]["head"]
            #    #dh_ub = h_i_res - JuMP.lower_bound(wm.var[:nw][n][:h][j])
            #    #dh_lb = h_i_res - JuMP.upper_bound(wm.var[:nw][n][:h][j])

            #    #mccormick_1[a][r] = JuMP.@constraint(wm.model, dh[a][r] >= dh_lb * x_r_a)
            #    #mccormick_2[a][r] = JuMP.@constraint(wm.model, dh[a][r] + h[j] >= dh_ub * x_r_a + h_i_res - dh_ub)
            #    #mccormick_3[a][r] = JuMP.@constraint(wm.model, dh[a][r] <= x_r_a * dh_ub)
            #    #mccormick_4[a][r] = JuMP.@constraint(wm.model, dh[a][r] + h[j] <= h_i_res + dh_lb * x_r_a - dh_lb)

            #    #diff_p_definition[a][r] = JuMP.@constraint(wm.model, diff_p_expr - inv(L) * dh[a][r] >= 0.0)
            #    #diff_n_definition[a][r] = JuMP.@constraint(wm.model, diff_n_expr + inv(L) * dh[a][r] >= 0.0)
            #elseif i in junction_ids
            #    #dh_ub = JuMP.upper_bound(wm.var[:nw][n][:h][i]) - JuMP.lower_bound(wm.var[:nw][n][:h][j])
            #    #dh_lb = JuMP.lower_bound(wm.var[:nw][n][:h][i]) - JuMP.upper_bound(wm.var[:nw][n][:h][j])

            #    #mccormick_1[a][r] = JuMP.@constraint(wm.model, dh[a][r] >= dh_lb * x_r_a)
            #    #mccormick_2[a][r] = JuMP.@constraint(wm.model, dh[a][r] - (h[i] - h[j]) >= dh_ub * x_r_a - dh_ub)
            #    #mccormick_3[a][r] = JuMP.@constraint(wm.model, dh[a][r] <= dh_ub * x_r_a)
            #    #mccormick_4[a][r] = JuMP.@constraint(wm.model, dh[a][r] - (h[i] - h[j]) <= dh_lb * x_r_a - dh_lb)

            #    #diff_p_definition[a][r] = JuMP.@constraint(wm.model, diff_p_expr - inv(L) * dh[a][r] >= 0.0)
            #    #diff_n_definition[a][r] = JuMP.@constraint(wm.model, diff_n_expr + inv(L) * dh[a][r] >= 0.0)
            #end

            # Add dual constraints of g variables.
            tmp_var_p = JuMP.@variable(wm.model, binary=true, start=0)
            JuMP.@constraint(wm.model, tmp_var_p <= x_r_a)
            JuMP.@constraint(wm.model, tmp_var_p <= x_dir)
            JuMP.@constraint(wm.model, tmp_var_p >= x_r_a + x_dir - 1)

            tmp_var_n = JuMP.@variable(wm.model, binary=true, start=0)
            JuMP.@constraint(wm.model, tmp_var_n <= x_r_a)
            JuMP.@constraint(wm.model, tmp_var_n <= 1 - x_dir)
            JuMP.@constraint(wm.model, tmp_var_n >= x_r_a + (1 - x_dir) - 1)

            #lambda_p_definition[a][r] = JuMP.@constraint(wm.model, sum(wm.var[:nw][n][:lambda_p][a][r]) == tmp_var_p)
            #lambda_p_definition[a][r] = JuMP.@constraint(wm.model, sum(wm.var[:nw][n][:lambda_n][a][r]) == tmp_var_n)
            lambda_p_definition[a][r] = JuMP.@constraint(wm.model, sum(wm.var[:nw][n][:lambda_p][a][r]) <= x_r_a)
            lambda_n_definition[a][r] = JuMP.@constraint(wm.model, sum(wm.var[:nw][n][:lambda_n][a][r]) <= x_r_a)
            lambda_p_definition[a][r] = JuMP.@constraint(wm.model, sum(wm.var[:nw][n][:lambda_p][a][r]) <= x_dir)
            lambda_p_definition[a][r] = JuMP.@constraint(wm.model, sum(wm.var[:nw][n][:lambda_n][a][r]) <= (1.0 - x_dir))
            lambda_p_definition[a][r] = JuMP.@constraint(wm.model, sum(wm.var[:nw][n][:lambda_p][a][r] .* qp_hat[a][r]) >= qp[a][r])
            lambda_n_definition[a][r] = JuMP.@constraint(wm.model, sum(wm.var[:nw][n][:lambda_n][a][r] .* qn_hat[a][r]) >= qn[a][r])

            wm.var[:nw][n][:fd] += (L * r_a) .* sum(wm.var[:nw][n][:lambda_p][a][r] .* gp_definition_rhs[a][r])
            wm.var[:nw][n][:fd] += (L * r_a) .* sum(wm.var[:nw][n][:lambda_n][a][r] .* gn_definition_rhs[a][r])
        end

        JuMP.@constraint(wm.model, sum(dhp_xr[a]) == dhp[a])
        JuMP.@constraint(wm.model, sum(dhn_xr[a]) == dhn[a])
    end

    wm.var[:nw][n][:fp] += sum(wm.var[:nw][n][:q_r][i] * res["head"] for (i, res) in ref(wm, n, :reservoirs))

    #for (i, reservoir) in wm.ref[:nw][n][:reservoirs]
    #    for (a, link) in filter(a -> i == a.second["node_fr"], wm.ref[:nw][n][:links])
    #        wm.var[:nw][n][:fp] += reservoir["head"] * (sum(qp[a]) - sum(qn[a]))
    #    end

    #    for (a, link) in filter(a -> i == a.second["node_to"], wm.ref[:nw][n][:links])
    #        wm.var[:nw][n][:fp] -= reservoir["head"] * (sum(qp[a]) - sum(qn[a]))
    #    end
    #end

    for (i, junction) in wm.ref[:nw][n][:junctions]
        wm.var[:nw][n][:fd] += wm.var[:nw][n][:h][i] * junction["demand"]
    end

    strong_duality = JuMP.@constraint(wm.model, wm.var[:nw][n][:fp] - wm.var[:nw][n][:fd] == 0.0)
end

"""
Implements a Benders-like algorithm that applies feasibility cuts to network
designs that are infeasible with respect to flow and/or head bounds.
"""
function ne_primal_dual(network_path::String, modifications_path::String,
                        mip_optimizer::JuMP.OptimizerFactory)
    # Read in the original network data.
    #alpha = network["options"]["headloss"] == "h-w" ? 1.852 : 2.0
    
    # Add the modifications to the network data and construct the MILPR model.
    network = WMs.parse_file(network_path)
    modifications = WMs.parse_file(modifications_path)
    InfrastructureModels.update_data!(network, modifications)

    # Build the relaxed network expansion problem.
    ne = WMs.build_generic_model(network, WMs.MILPRWaterModel, WMs.post_ne)
    add_system!(ne) # Add the primal-dual feasibility system.
    ne_solution = WMs.solve_generic_model(ne, mip_optimizer)

    ## EXTRA STUFF BELOW TO CHECK WHY SOLUTIONS ARE NOT INFEASIBLE...
    ## Get the corresponding (unbounded) CNLP solution.
    #ipopt = JuMP.with_optimizer(Ipopt.Optimizer, tol=1.0e-9, acceptable_tol=1.0e-9, print_level=0, linear_solver="mumps")
    #q, h, resistance_indices, cnlp_obj = get_cnlp_solution(ne, ipopt)

    #solve_without_bounds(ne, mip_optimizer)

    #println("FLOW SOLUTIONS")
    #for (a, link) in ne.ref[:nw][ne.cnw][:links]
    #    println(a, " ", JuMP.lower_bound(ne.var[:nw][ne.cnw][:q][a]), " ",
    #            JuMP.upper_bound(ne.var[:nw][ne.cnw][:q][a]), " ",
    #            q[a], " ", JuMP.value(ne.var[:nw][ne.cnw][:q][a]))
    #end

    ##println("")

    #println("LAMBDA SOLUTIONS")
    #for (a, link) in ne.ref[:nw][ne.cnw][:links]
    #    for r in 1:length(ne.ref[:nw][ne.cnw][:resistance][a])
    #        x_sel = r == resistance_indices[a] ? 1.0 : 0.0

    #        if x_sel == 1.0
    #            r_a = ne.ref[:nw][ne.cnw][:resistance][a][r]
    #            println(a, " + ", JuMP.value.(ne.var[:nw][ne.cnw][:lambda_p][a][r]), " ", sum(JuMP.value.(ne.var[:nw][ne.cnw][:lambda_p][a][r])))
    #            println(a, " - ", JuMP.value.(ne.var[:nw][ne.cnw][:lambda_n][a][r]), " ", sum(JuMP.value.(ne.var[:nw][ne.cnw][:lambda_n][a][r])))
    #        end
    #    end
    #end

    ##println("")

    #for (i, reservoir) in ne.ref[:nw][ne.cnw][:reservoirs]
    #    println(i, " ", reservoir["head"], " ", reservoir["head"])
    #end

    #println("HEAD SOLUTIONS FROM NE")
    #for (i, junction) in ne.ref[:nw][ne.cnw][:junctions]
    #    println(i, " ", JuMP.lower_bound(ne.var[:nw][ne.cnw][:h][i]), " ",
    #            JuMP.upper_bound(ne.var[:nw][ne.cnw][:h][i]), " ",
    #            h[i], " ", JuMP.value(ne.var[:nw][ne.cnw][:h][i]))
    #end

    #println("")

    #println("OBJECTIVE VALUES FROM NE")
    #println(cnlp_obj, " ", JuMP.value(ne.var[:nw][ne.cnw][:fp]), " ", JuMP.value(ne.var[:nw][ne.cnw][:fd]))
end

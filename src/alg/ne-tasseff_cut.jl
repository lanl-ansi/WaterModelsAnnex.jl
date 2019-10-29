import Dualization

function add_fixed_primal_model!(wm::AbstractWaterModel, model::JuMP.Model, n::Int=wm.cnw)
    num_outer_cuts = 25
    alpha = ref(wm, n, :alpha)

    link_ids = collect(WMs.ids(wm, n, :link))
    junction_ids = collect(WMs.ids(wm, n, :junction))
    reservoir_ids = collect(WMs.ids(wm, n, :reservoir))
    node_ids = vcat(junction_ids, reservoir_ids)
    resistances = ref(wm, n, :resistance)
    resistance_indices = Dict{Int, Int}(a => 0 for a in ids(wm, n, :link_ne))
    demand_sum = sum(ref(wm, n, :junction, k, "demand") for k in junction_ids)

    for (a, link) in wm.ref[:nw][n][:link_ne]
        x_res, r_id = findmax(JuMP.value.(wm.var[:nw][n][:x_res][a]))
        resistance_indices[a] = r_id
    end

    fp = JuMP.AffExpr(0.0)
    qp = JuMP.@variable(model, [a in link_ids], base_name="qp")
    qn = JuMP.@variable(model, [a in link_ids], base_name="qn")
    gp = JuMP.@variable(model, [a in link_ids], base_name="gp")
    gn = JuMP.@variable(model, [a in link_ids], base_name="gn")
    qr = JuMP.@variable(model, [i in reservoir_ids], base_name = "qr", lower_bound=0.0)

    for (a, link) in wm.ref[:nw][n][:link]
        L = link["length"]
        i = link["node_fr"]
        j = link["node_to"]
        r = resistance_indices[a]
        r_a = wm.ref[:nw][n][:resistance][a][r]

        qp_lower_bound = JuMP.@constraint(model, qp[a] >= 0.0)
        JuMP.set_name(qp_lower_bound, "qp_lower_bound[$(a)][$(r)]")
        qp_upper_bound = JuMP.@constraint(model, qp[a] <= demand_sum)
        JuMP.set_name(qp_upper_bound, "qp_upper_bound[$(a)][$(r)]")
        qp_hat = range(0.0, stop=demand_sum, length=num_outer_cuts)

        qn_lower_bound = JuMP.@constraint(model, qn[a] >= 0.0)
        JuMP.set_name(qn_lower_bound, "qn_lower_bound[$(a)][$(r)]")
        qn_upper_bound = JuMP.@constraint(model, qn[a] <= demand_sum)
        JuMP.set_name(qn_upper_bound, "qn_upper_bound[$(a)][$(r)]")
        qn_hat = range(0.0, stop=demand_sum, length=num_outer_cuts)

        # Add positive flow outer-approximation constraints.
        for k in 1:num_outer_cuts
            lhs_p = inv(1.0 + alpha) * qp_hat[k]^(1.0 + alpha) + qp_hat[k]^alpha * (qp[a] - qp_hat[k])
            gp_definition = JuMP.@constraint(model, lhs_p <= inv(L) * gp[a])
            JuMP.set_name(gp_definition, "gp_definition[$(a)][$(k)]")
        end

        # Add negative flow outer-approximation constraints.
        for k in 1:num_outer_cuts
            lhs_n = inv(1.0 + alpha) * qn_hat[k]^(1.0 + alpha) + qn_hat[k]^alpha * (qn[a] - qn_hat[k])
            gn_definition = JuMP.@constraint(model, lhs_n <= inv(L) * gn[a])
            JuMP.set_name(gn_definition, "gn_definition[$(a)][$(k)]")
        end

        fp += r_a * (gp[a] + gn[a])
    end

    fp -= sum(qr[i] * res["head"] for (i, res) in ref(wm, n, :reservoir))

    # Define flow conservation constraints.
    for (i, node) in wm.ref[:nw][n][:node]
        junc = ref(wm, n, :node_junction, i)
        arcs_fr = ref(wm, n, :node_arc_fr, i)
        arcs_to = ref(wm, n, :node_arc_to, i)
        res = ref(wm, n, :node_reservoir, i)
        tank = ref(wm, n, :node_tank, i)
        demands = Dict(k => ref(wm, n, :junction, k, "demand") for k in junc)

        c = JuMP.@constraint(model,
            sum(-qp[l] + qn[l] for (l, f, t) in arcs_fr) +
            sum( qp[l] - qn[l] for (l, f, t) in arcs_to) ==
            sum(-qr[rid] for rid in res) +
            sum(demand for (jid, demand) in demands))

        JuMP.set_name(c, "flow_conservation[$(i)]")
    end

    JuMP.@objective(model, MOI.MIN_SENSE, fp)
end

function add_primal_model!(wm::AbstractWaterModel, model::JuMP.Model, n::Int=wm.cnw; primal_dual::Bool=false)
    num_outer_cuts = 25
    alpha = wm.ref[:nw][n][:alpha]
    resistance_indices = Dict{Int, Int}(a => 0 for a in ids(wm, n, :link_ne))

    for (a, link) in wm.ref[:nw][n][:link_ne]
        x_res, r_id = findmax(JuMP.value.(wm.var[:nw][n][:x_res][a]))
        resistance_indices[a] = r_id
    end

    link_ids = collect(WMs.ids(wm, n, :link))
    junction_ids = collect(WMs.ids(wm, n, :junction))
    reservoir_ids = collect(WMs.ids(wm, n, :reservoir))
    node_ids = vcat(junction_ids, reservoir_ids)
    resistances = ref(wm, n, :resistance)

    fp = JuMP.AffExpr(0.0)
    qr = JuMP.@variable(model, [i in reservoir_ids], base_name = "qr")
    qp = JuMP.@variable(model, [a in link_ids, r in 1:length(resistances[a])], base_name="qp")
    qn = JuMP.@variable(model, [a in link_ids, r in 1:length(resistances[a])], base_name="qn")
    gp = JuMP.@variable(model, [a in link_ids, r in 1:length(resistances[a])], base_name="gp")
    gn = JuMP.@variable(model, [a in link_ids, r in 1:length(resistances[a])], base_name="gn")
    xr = JuMP.@variable(model, [a in link_ids, r in 1:length(resistances[a])], base_name="xr")

    for (a, link) in wm.ref[:nw][n][:link]
        L = link["length"]
        i = link["node_fr"]
        j = link["node_to"]

        for r in 1:length(resistances[a])
            r_a = wm.ref[:nw][n][:resistance][a][r]

            qp_a_r_ub = JuMP.upper_bound(wm.var[:nw][n][:qp_ne][a][r])
            qp_lower_bound = JuMP.@constraint(model, qp[a, r] >= 0.0)
            JuMP.set_name(qp_lower_bound, "qp_lower_bound[$(a)][$(r)]")

            qn_a_r_ub = JuMP.upper_bound(wm.var[:nw][n][:qn_ne][a][r])
            qn_lower_bound = JuMP.@constraint(model, qn[a, r] >= 0.0)
            JuMP.set_name(qn_lower_bound, "qn_lower_bound[$(a)][$(r)]")

            #xr_lower_bound = JuMP.@constraint(model, xr[a, r] >= 0.0)
            #xr_upper_bound = JuMP.@constraint(model, xr[a, r] <= 1.0)
            #JuMP.set_name(xr_lower_bound, "xr_lower_bound[$(a)][$(r)]")
            #JuMP.set_name(xr_upper_bound, "xr_upper_bound[$(a)][$(r)]")

            qp_upper_bound = JuMP.@constraint(model, qp[a, r] <= xr[a, r] * qp_a_r_ub)
            JuMP.set_name(qp_upper_bound, "qp_upper_bound[$(a)][$(r)]")
            qn_upper_bound = JuMP.@constraint(model, qn[a, r] <= xr[a, r] * qn_a_r_ub)
            JuMP.set_name(qn_upper_bound, "qn_upper_bound[$(a)][$(r)]")

            if r == resistance_indices[a]
                JuMP.fix(xr[a, r], 1.0, force=true)
                #JuMP.@constraint(model, xr[a, r] == 1.0)
            else
                JuMP.fix(xr[a, r], 0.0, force=true)
                #JuMP.@constraint(model, xr[a, r] == 0.0)
            end

            qp_hat = range(0.0, stop=qp_a_r_ub, length=num_outer_cuts)
            qn_hat = range(0.0, stop=qn_a_r_ub, length=num_outer_cuts)

            # Add positive flow outer-approximation constraints.
            for k in 1:num_outer_cuts
                lhs_p = inv(1.0 + alpha) * qp_hat[k]^(1.0 + alpha) + qp_hat[k]^alpha * (qp[a, r] - qp_hat[k])
                gp_definition = JuMP.@constraint(model, lhs_p <= inv(L) * gp[a, r])
                JuMP.set_name(gp_definition, "gp_definition[$(a)][$(r)][$(k)]")
            end

            # Add negative flow outer-approximation constraints.
            for k in 1:num_outer_cuts
                lhs_n = inv(1.0 + alpha) * qn_hat[k]^(1.0 + alpha) + qn_hat[k]^alpha * (qn[a, r] - qn_hat[k])
                gn_definition = JuMP.@constraint(model, lhs_n <= inv(L) * gn[a, r])
                JuMP.set_name(gn_definition, "gn_definition[$(a)][$(r)][$(k)]")
            end

            fp += r_a * (gp[a, r] + gn[a, r])
        end
    end

    # Define flow conservation constraints.
    fp -= sum(qr[i] * res["head"] for (i, res) in ref(wm, n, :reservoir))

    for (i, node) in wm.ref[:nw][n][:node]
        junc = ref(wm, n, :node_junction, i)
        arcs_fr = ref(wm, n, :node_arc_fr, i)
        arcs_to = ref(wm, n, :node_arc_to, i)
        res = ref(wm, n, :node_reservoir, i)
        tank = ref(wm, n, :node_tank, i)
        demands = Dict(k => ref(wm, n, :junction, k, "demand") for k in junc)

        c = JuMP.@constraint(model,
            sum(-sum(qp[l, r] - qn[l, r] for r in 1:length(resistances[l])) for (l, f, t) in arcs_fr) +
            sum( sum(qp[l, r] - qn[l, r] for r in 1:length(resistances[l])) for (l, f, t) in arcs_to) ==
            sum(-qr[rid] for rid in res) +
            sum(demand for (jid, demand) in demands))

        JuMP.set_name(c, "flow_conservation[$(i)]")
    end

    #delta = JuMP.@variable(model, base_name="Delta", lower_bound=0.0)
    JuMP.@constraint(model, fp - JuMP.objective_function(model) == 0.0)
    #JuMP.@objective(model, MOI.MIN_SENSE, delta)
end

function add_dual_bounds!(wm::AbstractWaterModel, model::JuMP.Model, n::Int=wm.cnw)
    for (i, node) in wm.ref[:nw][n][:node]
        var = JuMP.variable_by_name(model, "dual_flow_conservation[$(i)]_1")
        lb = JuMP.lower_bound(WMs.var(wm, n, :h, i))
        JuMP.@constraint(model, -var >= lb)
        ub = JuMP.upper_bound(WMs.var(wm, n, :h, i))
        JuMP.@constraint(model, -var <= ub)
    end
end

function get_dual_model(model::JuMP.Model)
    dual_model = Dualization.dualize(model, dual_names=Dualization.DualNames("dual_", ""))
end

function get_infeasibility_ray(wm::AbstractWaterModel, model::JuMP.Model, n::Int=wm.cnw)
    c_leq = JuMP.all_constraints(model, JuMP.AffExpr, MOI.LessThan{Float64})
    c_geq = JuMP.all_constraints(model, JuMP.AffExpr, MOI.GreaterThan{Float64})
    c_eq = JuMP.all_constraints(model, JuMP.AffExpr, MOI.EqualTo{Float64})
    v_leq = JuMP.all_constraints(model, JuMP.VariableRef, MOI.LessThan{Float64})
    v_geq = JuMP.all_constraints(model, JuMP.VariableRef, MOI.GreaterThan{Float64})
    v_eq = JuMP.all_constraints(model, JuMP.VariableRef, MOI.EqualTo{Float64})

    lhs = JuMP.AffExpr(0.0)

    for c in c_geq
        con = JuMP.constraint_object(c)
        constant = con.set.lower
        lhs += JuMP.dual(c) * constant
    end

    for c in c_eq
        con = JuMP.constraint_object(c)
        constant = con.set.value
        vars = collect(keys(con.func.terms))
        x_vars = filter(x -> occursin("xr", JuMP.name(x)), vars)

        if length(x_vars) > 0
            var_name = JuMP.name(x_vars[1])
            indices = split(split(var_name, "[")[2], "]")[1]
            a = parse(Int, split(indices, ",")[1])
            r = parse(Int, split(indices, ",")[2])
            x_r_a = wm.var[:nw][n][:x_res][a][r]
            lhs -= JuMP.dual(c) * x_r_a * con.func.terms[x_vars[1]]
        else
            lhs += JuMP.dual(c) * constant
        end
    end

    for c in v_leq
        con = JuMP.constraint_object(c)
        constant = con.set.upper
        lhs += JuMP.dual(c) * constant
    end

    for c in v_geq
        con = JuMP.constraint_object(c)
        constant = con.set.lower
        lhs += JuMP.dual(c) * constant
    end

    for c in v_eq
        con = JuMP.constraint_object(c)
        constant = con.set.value
        lhs += JuMP.dual(c) * constant
    end

    for c in c_leq
        con = JuMP.constraint_object(c)
        constant = con.set.upper

        vars = collect(keys(con.func.terms))
        x_vars = filter(x -> occursin("xr", JuMP.name(x)), vars)

        if length(x_vars) > 0
            var_name = JuMP.name(x_vars[1])
            indices = split(split(var_name, "[")[2], "]")[1]
            a = parse(Int, split(indices, ",")[1])
            r = parse(Int, split(indices, ",")[2])
            x_r_a = wm.var[:nw][n][:x_res][a][r]
            lhs -= JuMP.dual(c) * x_r_a * con.func.terms[x_vars[1]]
        else
            lhs += JuMP.dual(c) * constant
        end
    end

    return lhs
end

function add_tasseff_cut!(wm::AbstractWaterModel, optimizer::JuMP.OptimizerFactory, n::Int=wm.cnw)
    model = JuMP.Model(JuMP.with_optimizer(optimizer))
    add_fixed_primal_model!(wm, model, n)
    model = get_dual_model(model)
    add_dual_bounds!(wm, model, n)
    add_primal_model!(wm, model, n)
    JuMP.optimize!(model, optimizer)

    # Get and add the infeasibility ray. 
    lhs = get_infeasibility_ray(wm, model, n)
    JuMP.@constraint(wm.model, lhs >= 0.0)

    #for (i, node) in wm.ref[:nw][n][:node]
    #    #con = JuMP.constraint_by_name(model, "flow_conservation[$(i)]")
    #    #println(i, " ", JuMP.shadow_price(con))

    #    var = JuMP.variable_by_name(model, "dual_flow_conservation[$(i)]_1")
    #    println(i, " ", JuMP.value(var))
    #    delta = JuMP.variable_by_name(model, "Delta[$(i)]")
    #end

    #add_primal_model!(wm, primal_model, n, primal_dual=false)
    #add_primal_model!(wm, model, n, primal_dual=true)

    ##get_flow_solution(wm, model, n)
    ##delta = JuMP.variable_by_name(model, "Delta")
    ##println(JuMP.value(delta))


    #alpha = wm.ref[:nw][n][:alpha]
 
    #resistance_indices = Dict{Int, Int}(a => 0 for a in ids(wm, n, :link_ne))
    #for (a, link) in wm.ref[:nw][n][:link_ne]
    #    x_res, r_id = findmax(JuMP.value.(wm.var[:nw][n][:x_res][a]))
    #    resistance_indices[a] = r_id
    #end

    #link_ids = collect(WMs.ids(wm, n, :link))
    #junction_ids = collect(WMs.ids(wm, n, :junction))
    #reservoir_ids = collect(WMs.ids(wm, n, :reservoir))
    #node_ids = vcat(junction_ids, reservoir_ids)
    #resistances = ref(wm, n, :resistance)

    #fp = JuMP.AffExpr(0.0)
    #fd = JuMP.AffExpr(0.0)

    #mccormick_1 = Dict{Int, JuMP.ConstraintRef}()
    #mccormick_2 = Dict{Int, JuMP.ConstraintRef}()
    #mccormick_3 = Dict{Int, JuMP.ConstraintRef}()
    #mccormick_4 = Dict{Int, JuMP.ConstraintRef}()

    #qp_hat = Dict{Int, Dict{Int, Array{Float64}}}(a => Dict{Int, Array{Float64}}(r => [0] for r in 1:length(resistances[a])) for a in link_ids)
    #qn_hat = Dict{Int, Dict{Int, Array{Float64}}}(a => Dict{Int, Array{Float64}}(r => [0] for r in 1:length(resistances[a])) for a in link_ids)

    #qp_upper_bound = Dict{Int, Dict{Int, JuMP.ConstraintRef}}(a => Dict{Int, JuMP.ConstraintRef}() for a in link_ids)
    #qp_lower_bound = Dict{Int, Dict{Int, JuMP.ConstraintRef}}(a => Dict{Int, JuMP.ConstraintRef}() for a in link_ids)
    #qn_upper_bound = Dict{Int, Dict{Int, JuMP.ConstraintRef}}(a => Dict{Int, JuMP.ConstraintRef}() for a in link_ids)
    #qn_lower_bound = Dict{Int, Dict{Int, JuMP.ConstraintRef}}(a => Dict{Int, JuMP.ConstraintRef}() for a in link_ids)
    #xr_upper_bound = Dict{Int, Dict{Int, JuMP.ConstraintRef}}(a => Dict{Int, JuMP.ConstraintRef}() for a in link_ids)
    #xr_lower_bound = Dict{Int, Dict{Int, JuMP.ConstraintRef}}(a => Dict{Int, JuMP.ConstraintRef}() for a in link_ids)
    #gp_definition = Dict{Int, Dict{Int, Dict{Int, JuMP.ConstraintRef}}}(a => Dict{Int, Dict{Int, JuMP.ConstraintRef}}(r => Dict{Int, JuMP.ConstraintRef}() for r in 1:length(resistances[a])) for a in link_ids)
    #gn_definition = Dict{Int, Dict{Int, Dict{Int, JuMP.ConstraintRef}}}(a => Dict{Int, Dict{Int, JuMP.ConstraintRef}}(r => Dict{Int, JuMP.ConstraintRef}() for r in 1:length(resistances[a])) for a in link_ids)

    ##gp_definition_rhs = Dict{Int, Array{Float64}}()
    ##gn_definition_rhs = Dict{Int, Array{Float64}}()

    #num_outer_cuts_p = Dict{Int, Dict{Int, Int}}(a => Dict{Int, Int}(r => 0 for r in 1:length(resistances[a])) for a in link_ids)
    #num_outer_cuts_n = Dict{Int, Dict{Int, Int}}(a => Dict{Int, Int}(r => 0 for r in 1:length(resistances[a])) for a in link_ids)

    ##lambda_p = Dict{Int, Array{JuMP.VariableRef}}()
    ##lambda_n = Dict{Int, Array{JuMP.VariableRef}}()

    #diff_p_definition = Dict{Int, JuMP.ConstraintRef}()
    #diff_n_definition = Dict{Int, JuMP.ConstraintRef}()

    ##lambda_p_definition = Dict{Int, JuMP.ConstraintRef}()
    ##lambda_n_definition = Dict{Int, JuMP.ConstraintRef}()

    #xr = JuMP.@variable(model, [a in link_ids, r in 1:length(resistances[a])], base_name="xr[$a]")
    #qp = JuMP.@variable(model, [a in link_ids, r in 1:length(resistances[a])], base_name="qp[$a]")
    #qn = JuMP.@variable(model, [a in link_ids, r in 1:length(resistances[a])], base_name="qn[$a]")
    #gp = JuMP.@variable(model, [a in link_ids, r in 1:length(resistances[a])], base_name="gp[$a]")
    #gn = JuMP.@variable(model, [a in link_ids, r in 1:length(resistances[a])], base_name="gn[$a]")
    #qr = JuMP.@variable(model, [i in reservoir_ids], base_name = "qr")

    ##dh = JuMP.@variable(model, [a in link_ids], base_name="dh[$a]")
    ##h = JuMP.@variable(model, [i in node_ids], base_name = "h")

    #for (a, link) in wm.ref[:nw][n][:link]
    #    L = link["length"]
    #    i = link["node_fr"]
    #    j = link["node_to"]

    #    for r in 1:length(resistances[a])
    #        r_a = wm.ref[:nw][n][:resistance][a][r]

    #        if r == resistance_indices[a]
    #            JuMP.fix(xr[a, r], 1.0, force=true)
    #        else
    #            JuMP.fix(xr[a, r], 0.0, force=true)
    #        end

    #        xr_lower_bound[a][r] = JuMP.@constraint(model, xr[a, r] >= 0.0)
    #        xr_upper_bound[a][r] = JuMP.@constraint(model, xr[a, r] <= 1.0)
    #        JuMP.set_name(xr_lower_bound[a][r], "xr_lower_bound[$(a)][$(r)]")
    #        JuMP.set_name(xr_upper_bound[a][r], "xr_upper_bound[$(a)][$(r)]")

    #        qp_a_r_ub = JuMP.upper_bound(wm.var[:nw][n][:qp_ne][a][r])
    #        qp_lower_bound[a][r] = JuMP.@constraint(model, qp[a, r] >= 0.0)
    #        qp_upper_bound[a][r] = JuMP.@constraint(model, qp[a, r] <= xr[a, r] * qp_a_r_ub)
    #        JuMP.set_name(qp_lower_bound[a][r], "qp_lower_bound[$(a)][$(r)]")
    #        JuMP.set_name(qp_upper_bound[a][r], "qp_upper_bound[$(a)][$(r)]")

    #        qn_a_r_ub = JuMP.upper_bound(wm.var[:nw][n][:qn_ne][a][r])
    #        qn_lower_bound[a][r] = JuMP.@constraint(model, qn[a, r] >= 0.0)
    #        qn_upper_bound[a][r] = JuMP.@constraint(model, qn[a, r] <= xr[a, r] * qn_a_r_ub)
    #        JuMP.set_name(qn_lower_bound[a][r], "qn_lower_bound[$(a)][$(r)]")
    #        JuMP.set_name(qn_upper_bound[a][r], "qn_upper_bound[$(a)][$(r)]")

    #        num_outer_cuts_p[a][r] = 25
    #        qp_hat[a][r] = range(0.0, stop=qp_a_r_ub, length=num_outer_cuts_p[a][r])
    #        num_outer_cuts_n[a][r] = 25
    #        qn_hat[a][r] = range(0.0, stop=qn_a_r_ub, length=num_outer_cuts_n[a][r])

    #        # Add positive flow outer-approximation constraints.
    #        for k in 1:num_outer_cuts_p[a][r]
    #            rhs_p = r_a * (1.0 - inv(1.0 + alpha)) * qp_hat[a][r][k]^(1.0 + alpha)
    #            lhs_p = r_a * qp_hat[a][r][k]^alpha * qp[a, r] - gp[a, r]
    #            gp_definition[a][r][k] = JuMP.@constraint(model, lhs_p <= rhs_p)
    #            JuMP.set_name(gp_definition[a][r][k], "gp_definition[$(a)][$(r)][$(k)]")
    #        end

    #        # Add negative flow outer-approximation constraints.
    #        for k in 1:num_outer_cuts_n[a][r]
    #            rhs_n = r_a * (1.0 - inv(1.0 + alpha)) * qn_hat[a][r][k]^(1.0 + alpha)
    #            lhs_n = r_a * qn_hat[a][r][k]^alpha * qn[a, r] - gn[a, r]
    #            gn_definition[a][r][k] = JuMP.@constraint(model, lhs_n <= rhs_n)
    #            JuMP.set_name(gn_definition[a][r][k], "gn_definition[$(a)][$(r)][$(k)]")
    #        end

    #        fp += L * (gp[a, r] + gn[a, r])
    #    end

    ##    gp_definition[a] = Dict{Int, JuMP.ConstraintRef}()
    ##    gn_definition[a] = Dict{Int, JuMP.ConstraintRef}()
    ##    gp_definition_rhs[a] = zeros(num_outer_cuts_p[a])
    ##    gn_definition_rhs[a] = zeros(num_outer_cuts_n[a])

    ##    # Add positive flow outer-approximation constraints.
    ##    for k in 1:num_outer_cuts_p[a]
    ##        rhs_p = (1.0 - inv(1.0 + alpha)) * qp_hat[a][k]^(1.0 + alpha)
    ##        lhs_p = qp_hat[a][k]^alpha * qp[a] - gp[a]
    ##        gp_definition[a][k] = JuMP.@constraint(model, lhs_p <= rhs_p)
    ##        gp_definition_rhs[a][k] = rhs_p
    ##        JuMP.set_name(gp_definition[a][k], "gp_definition[$(a)][$(k)]")
    ##    end

    ##    # Add negative flow outer-approximation constraints.
    ##    for k in 1:num_outer_cuts_n[a]
    ##        rhs_n = (1.0 - inv(1.0 + alpha)) * qn_hat[a][k]^(1.0 + alpha)
    ##        lhs_n = qn_hat[a][k]^alpha * qn[a] - gn[a]
    ##        gn_definition[a][k] = JuMP.@constraint(model, lhs_n <= rhs_n)
    ##        gn_definition_rhs[a][k] = rhs_n
    ##        JuMP.set_name(gn_definition[a][k], "gn_definition[$(a)][$(k)]")
    ##    end

    ##    #lambda_p[a] = JuMP.@variable(model, [k in 1:num_outer_cuts_p[a]], lower_bound=0.0, base_name="lambda_p[$(a)]")
    ##    #lambda_n[a] = JuMP.@variable(model, [k in 1:num_outer_cuts_n[a]], lower_bound=0.0, base_name="lambda_n[$(a)]")

    ##    ## Add dual constraints of flow variables.
    ##    #diff_p_expr = r_a * sum(qp_hat[a].^alpha .* lambda_p[a])
    ##    #diff_n_expr = r_a * sum(qn_hat[a].^alpha .* lambda_n[a])

    ##    #dh_ub = JuMP.upper_bound(wm.var[:nw][n][:h][i]) - JuMP.lower_bound(wm.var[:nw][n][:h][j])
    ##    #dh_lb = JuMP.lower_bound(wm.var[:nw][n][:h][i]) - JuMP.upper_bound(wm.var[:nw][n][:h][j])

    ##    #mccormick_1[a] = JuMP.@constraint(model, dh[a] >= dh_lb)
    ##    #mccormick_2[a] = JuMP.@constraint(model, dh[a] - (h[i] - h[j]) >= 0.0)
    ##    #mccormick_3[a] = JuMP.@constraint(model, dh[a] <= dh_ub)
    ##    #mccormick_4[a] = JuMP.@constraint(model, dh[a] - (h[i] - h[j]) <= 0.0)

    ##    #diff_p_definition[a] = JuMP.@constraint(model, diff_p_expr - inv(L) * dh[a] >= 0.0)
    ##    #diff_n_definition[a] = JuMP.@constraint(model, diff_n_expr + inv(L) * dh[a] >= 0.0)

    ##    ## Add dual constraints of g variables.
    ##    #lambda_p_definition[a] = JuMP.@constraint(model, sum(lambda_p[a]) <= 1.0)
    ##    #lambda_n_definition[a] = JuMP.@constraint(model, sum(lambda_n[a]) <= 1.0)

    ##    #fd += L * r_a * sum(lambda_p[a] .* gp_definition_rhs[a])
    ##    #fd += L * r_a * sum(lambda_n[a] .* gn_definition_rhs[a])
    #end

    ## Define flow conservation constraints.
    #fp -= sum(qr[i] * res["head"] for (i, res) in ref(wm, n, :reservoir))
    #flow_conservation = Dict{Int, JuMP.ConstraintRef}()

    #for (i, node) in wm.ref[:nw][n][:node]
    #    junc = ref(wm, n, :node_junction, i)
    #    arcs_fr = ref(wm, n, :node_arc_fr, i)
    #    arcs_to = ref(wm, n, :node_arc_to, i)
    #    res = ref(wm, n, :node_reservoir, i)
    #    tank = ref(wm, n, :node_tank, i)
    #    demands = Dict(k => ref(wm, n, :junction, k, "demand") for k in junc)

    #    flow_conservation[i] = JuMP.@constraint(model,
    #        sum(sum(-qp[l, r] + qn[l, r] for r in 1:length(resistances[l])) for (l, f, t) in arcs_fr) +
    #        sum(sum( qp[l, r] - qn[l, r] for r in 1:length(resistances[l])) for (l, f, t) in arcs_to) ==
    #        sum(-qr[rid] for rid in res) +
    #        sum(demand for (jid, demand) in demands))

    #    JuMP.set_name(flow_conservation[i], "flow_conservation[$(i)]")
    #end

    #h_lower_bound = Dict{Int, JuMP.ConstraintRef}()
    #h_upper_bound = Dict{Int, JuMP.ConstraintRef}()

    ###for (i, node) in wm.ref[:nw][n][:node]
    ###    h_i_lb = JuMP.lower_bound(wm.var[:nw][n][:h][i])
    ###    h_i_ub = JuMP.upper_bound(wm.var[:nw][n][:h][i])
    ###    h_lower_bound[i] = JuMP.@constraint(model, h[i] >= h_i_lb)
    ###    h_upper_bound[i] = JuMP.@constraint(model, h[i] <= h_i_ub)
    ###end

    ###for (i, junction) in wm.ref[:nw][n][:junction]
    ###    fd += h[i] * junction["demand"]
    ###end

    #JuMP.@objective(model, MOI.MIN_SENSE, fp)

    ###strong_duality = JuMP.@constraint(model, fp - fd == 0.0)

    ##println(JuMP.dual(flow_conservation[3]))
    ##println(JuMP.objective_value(model))
    #
    #dual_model = Dualization.dualize(model, dual_names=Dualization.DualNames("dual_", ""))

    #for (i, node) in wm.ref[:nw][n][:node]
    #    var = JuMP.variable_by_name(dual_model, "dual_flow_conservation[$(i)]_1")
    #    lb = JuMP.lower_bound(WMs.var(wm, n, :h, i))
    #    ub = JuMP.upper_bound(WMs.var(wm, n, :h, i))
    #    h_lower_bound[i] = JuMP.@constraint(dual_model, var >= lb)
    #    h_upper_bound[i] = JuMP.@constraint(dual_model, var <= ub)
    #end

    #fd = JuMP.objective_function(dual_model)
    #JuMP.@constraint(dual_model, fp == fd)
    #JuMP.optimize!(dual_model, optimizer)

    #for (i, node) in wm.ref[:nw][n][:node]
    #    var = JuMP.variable_by_name(dual_model, "dual_flow_conservation[$(i)]_1")
    #    println(i, " ", JuMP.value(var))
    #end

    ##println(dual_model)

    ##println(JuMP.value.(qp))

    ###if JuMP.dual_status(model) == MOI.INFEASIBILITY_CERTIFICATE
    ###    gt_cut_lhs = JuMP.AffExpr(0.0)

    ###    for (a, link) in wm.ref[:nw][n][:link_ne]
    ###        L = link["length"]
    ###        i = link["node_fr"]
    ###        j = link["node_to"]

    ###        r = resistance_indices[a]
    ###        r_a = wm.ref[:nw][n][:resistance][a][r]
    ###        x_r_a = wm.var[:nw][n][:x_res][a][r]

    ###        q_p_a_ub = JuMP.upper_bound(wm.var[:nw][n][:qp_ne][a][r])
    ###        q_n_a_ub = JuMP.upper_bound(wm.var[:nw][n][:qn_ne][a][r])

    ###        dh_ub = JuMP.upper_bound(wm.var[:nw][n][:h][i]) - JuMP.lower_bound(wm.var[:nw][n][:h][j])
    ###        dh_lb = JuMP.lower_bound(wm.var[:nw][n][:h][i]) - JuMP.upper_bound(wm.var[:nw][n][:h][j])

    ###        gt_cut_lhs += JuMP.dual(mccormick_1[a]) * (dh_lb * x_r_a)
    ###        gt_cut_lhs += JuMP.dual(mccormick_2[a]) * (dh_ub * x_r_a - dh_ub)
    ###        gt_cut_lhs += JuMP.dual(mccormick_3[a]) * (dh_ub * x_r_a)
    ###        gt_cut_lhs += JuMP.dual(mccormick_4[a]) * (dh_lb * x_r_a - dh_lb)
    ###        gt_cut_lhs += JuMP.dual(lambda_p_definition[a]) * x_r_a
    ###        gt_cut_lhs += JuMP.dual(lambda_n_definition[a]) * x_r_a

    ###        for k in 1:num_outer_cuts_p[a]
    ###            gt_cut_lhs += JuMP.dual(gp_definition[a][k]) * gp_definition_rhs[a][k]
    ###        end

    ###        for k in 1:num_outer_cuts_n[a]
    ###            gt_cut_lhs += JuMP.dual(gn_definition[a][k]) * gn_definition_rhs[a][k]
    ###        end

    ###        gt_cut_lhs += JuMP.dual(qp_upper_bound[a]) * q_p_a_ub * x_r_a
    ###        gt_cut_lhs += JuMP.dual(qn_upper_bound[a]) * q_n_a_ub * x_r_a
    ###    end

    ###    for (i, nodes) in wm.ref[:nw][n][:node]
    ###        gt_cut_lhs += JuMP.dual(h_lower_bound[i]) * JuMP.lower_bound(wm.var[:nw][n][:h][i])
    ###        gt_cut_lhs += JuMP.dual(h_upper_bound[i]) * JuMP.upper_bound(wm.var[:nw][n][:h][i])

    ###        if i in keys(wm.ref[:nw][n][:junction])
    ###            junction = wm.ref[:nw][n][:junction][i]
    ###            gt_cut_lhs += JuMP.dual(flow_conservation[i]) * junction["demand"]
    ###        end
    ###    end

    ###    #add_no_good_cut!(wm, optimizer, n)
    ###    c = JuMP.@constraint(wm.model, gt_cut_lhs >= 0.0)
    ###    #println(c)
    ###else
    ###    #add_no_good_cut!(wm, optimizer, n)
    ###end
 
    #add_no_good_cut!(wm, optimizer, n)
end

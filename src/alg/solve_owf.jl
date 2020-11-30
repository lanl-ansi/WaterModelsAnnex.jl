function solve_obbt(network_path::String, obbt_optimizer; time_limit::Float64 = 3600.0)
    # Read in the original network data.
    network = WM.parse_file(network_path)

    # Tighten bounds of variables in the network.
    WM.run_obbt_owf!(
        network, obbt_optimizer; model_type = LRDWaterModel,
        solve_relaxed = false, time_limit = time_limit,
        ext = Dict(:pipe_breakpoints => 10, :pump_breakpoints => 10))

    # Get tightened network data.
    return network
end


function solve_owf(network_path::String, obbt_optimizer, owf_optimizer, nlp_optimizer)
    # Tighten the bounds in the network.
    network = solve_obbt(network_path, obbt_optimizer)
    result = solve_owf(network_path, network, obbt_optimizer, owf_optimizer, nlp_optimizer)
    return result
end


function compute_pairwise_cuts(network::Dict{String, Any}, obbt_optimizer)
    # Get pairwise cutting planes from the network-relaxed problem.
    ext = Dict{Symbol, Any}(:pipe_breakpoints => 10, :pump_breakpoints => 10)
    wm = instantiate_model(network, LRDWaterModel, build_owf; ext = ext)
    WM.JuMP.set_optimizer(wm.model, obbt_optimizer)
    problem_sets = WM._get_pairwise_problem_sets(wm; nw = wm.cnw)
    return WM._compute_pairwise_cuts!(wm, problem_sets)
end


function construct_owf_model(network::Dict{String, Any}, owf_optimizer)
    # Construct the OWF model.
    network_mn = WM.make_multinetwork(network)
    ext = Dict(:pipe_breakpoints => 5, :pump_breakpoints => 5)
    wm = WM.instantiate_model(network_mn, LRDWaterModel, build_mn_owf; ext = ext)

    # Constrain an auxiliary objective variable by the objective function.
    objective_function = WM.JuMP.objective_function(wm.model)
    objective_var = WM.JuMP.@variable(wm.model, base_name = "obj_aux", lower_bound = 0.0)
    WM.JuMP.@objective(wm.model, WM._MOI.MIN_SENSE, objective_var)
    WM.JuMP.@constraint(wm.model, objective_function <= objective_var)

    # Set the optimizer and other important solver parameters.
    WM.JuMP.set_optimizer(wm.model, owf_optimizer)
    WM._MOI.set(wm.model, WM._MOI.NumberOfThreads(), 1)

    # Set the pump statuses for the Van Zyl network.
    #_set_van_zyl_statuses(wm)

    # Return the WaterModels object.
    return wm
end


function construct_relaxed_owf_model(network::Dict{String, Any}, nlp_optimizer)
    # Construct the relaxed OWF model.
    network_mn = WM.make_multinetwork(network)
    ext = Dict(:pipe_breakpoints => 5, :pump_breakpoints => 5)
    wm = WM.instantiate_model(network_mn, LRDWaterModel, build_mn_owf; ext = ext)
    WM.relax_all_binary_variables!(wm)

    # Set the optimizer and other important solver parameters.
    WM.JuMP.set_optimizer(wm.model, nlp_optimizer)

    # Return the WaterModels object.
    return wm
end


function add_pairwise_cuts(wm::AbstractWaterModel, cuts::Array{WM._PairwiseCut, 1})
    # Add the pairwise cuts obtained from the relaxed problem to the OWF problem.
    for nw_id in WM.nw_ids(wm)
        # Use the same cuts for all subnetworks of the multinetwork.
        map(x -> x.variable_index_1.network_index = nw_id, cuts)
        map(x -> x.variable_index_2.network_index = nw_id, cuts)

        # Add the collection of pairwise cuts for the subnetwork.
        WM._add_pairwise_cuts!(wm, cuts)
    end
end


function _get_lazy_cut_callback(wm::AbstractWaterModel)
    return function lazy_cut_callback(cb_data)
    end
end


function solve_owf(network_path::String, network, obbt_optimizer, owf_optimizer, nlp_optimizer)
    wm = construct_owf_model(network, owf_optimizer)
    pairwise_cuts = compute_pairwise_cuts(network, obbt_optimizer)
    add_pairwise_cuts(wm, pairwise_cuts)

    #lazy_cut_callback = _get_lazy_cut_callback(wm)
    wntr_network = WMA.initialize_wntr_network(wm.data)

    # Solve the convex, continuously-relaxed optimal water flow problem.
    wm_nlp = construct_relaxed_owf_model(network, nlp_optimizer)
    nlp_result = WM.optimize_model!(wm_nlp)

    num_infeasible_solutions = 0 # Number of integer *infeasible* solutions.
    objective_comparison_table = [] # Compares relaxed versus true objectives.
    num_simulations_executed = 0

    function lazy_cut_callback(cb_data) # Define the lazy cut callback function.
        # Populate the solution of wm_nlp to use in the feasibility check.
        _populate_solution!(cb_data, wm, wm_nlp)
        WMA.update_wntr_controls(wntr_network, wm.data, wm_nlp.solution, Float64(wm.data["time_step"]))
        wntr_result = WMA.simulate_wntr(wntr_network)
        num_simulations_executed += 1

        # Get infeasibilities associated with tanks
        tank_infeasibilities = _calc_tank_infeasibilities(wm, wntr_result)

        if sum(sum(abs.(array)) for array in values(tank_infeasibilities)) > 0.0
            # Collect the current integer solution into "zero" and "one" buckets.
            vars = _get_pump_indicators(wm) # All relevant component status variables.
            zero_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 0.0, vars)

            # If the solution is not feasible (according to a comparison with WNTR), add a no-good cut.
            con = WM.JuMP.@build_constraint(sum(zero_vars) >= 1.0)
            WM._MOI.submit(wm.model, WM._MOI.LazyConstraint(cb_data), con)

            # Increase the number of infeasible solutions encountered by one.
            num_infeasible_solutions += 1
        else
            # Get the first network index corresponding to infeasibility.
            min_nw = calc_first_nw_infeasibility(wm, wntr_result)

            if min_nw != nothing # There is an infeasibility somewhere.
                # If there is an infeasibility, add a no-good cut up to that network index.
                vars = _get_indicator_variables_to_nw(wm, min_nw)
                zero_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 0.0, vars)
                one_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 1.0, vars)

                # If the solution is not feasible (according to a comparison with WNTR), add a no-good cut.
                con = WM.JuMP.@build_constraint(sum(zero_vars) - sum(one_vars) >= 1.0 - length(one_vars))
                WM._MOI.submit(wm.model, WM._MOI.LazyConstraint(cb_data), con)
                num_infeasible_solutions += 1
            else # There are no infeasibilities.
                # Collect the current integer solution into "zero" and "one" buckets.
                vars = _get_indicator_variables(wm) # All relevant component status variables.
                zero_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 0.0, vars)
                one_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 1.0, vars)

                # Compute the objectives of the relaxation and WNTR solution.
                objective_var = WM.JuMP.variable_by_name(wm.model, "obj_aux")
                relaxed_objective = WM.JuMP.callback_value(cb_data, objective_var)
                true_objective = _calc_wntr_objective(wm, wntr_result)

                # Add a cut to raise the objective for the solution to the true objective.
                bin_expr = true_objective * (length(one_vars) - sum(one_vars) + sum(zero_vars))
                con = WM.JuMP.@build_constraint(objective_var >= true_objective - bin_expr)
                WM._MOI.submit(wm.model, WM._MOI.LazyConstraint(cb_data), con)
                push!(objective_comparison_table, (relaxed_objective, true_objective))
            end
        end
    end

    # Register the lazy cut callback with the JuMP modeling object.
    WM._MOI.set(wm.model, WM._MOI.LazyConstraintCallback(), lazy_cut_callback)

    #function user_cut_callback(cb_data) # Define the user cut callback function.
    #    for nw in WM.nw_ids(wm)
    #        exponent = WM.ref(wm, nw, :alpha)
    #        pipe_ids = collect(WM.ids(wm, nw, :pipe))
    #        L = [WM.ref(wm, nw, :pipe, i)["length"] for i in pipe_ids]
    #        r = [WM.ref(wm, nw, :resistance, i)[1] for i in pipe_ids]

    #        qp, qn = WM.var(wm, nw, :qp_pipe), WM.var(wm, nw, :qn_pipe)
    #        qp_vals = [max(0.0, WM.JuMP.callback_value(cb_data, qp[i])) for i in pipe_ids]
    #        qn_vals = [max(0.0, WM.JuMP.callback_value(cb_data, qn[i])) for i in pipe_ids]

    #        dhp, dhn = WM.var(wm, nw, :dhp_pipe), WM.var(wm, nw, :dhn_pipe)
    #        dhp_vals = [max(0.0, WM.JuMP.callback_value(cb_data, dhp[i])) for i in pipe_ids]
    #        dhn_vals = [max(0.0, WM.JuMP.callback_value(cb_data, dhn[i])) for i in pipe_ids]

    #        dhp_ests = L .* r .* qp_vals.^(exponent)
    #        dhn_ests = L .* r .* qn_vals.^(exponent)

    #        for i in 1:length(pipe_ids)
    #            if abs(dhp_vals[i] - dhp_ests[i]) > 1.0e-1
    #                y = WM.var(wm, nw, :y_pipe, pipe_ids[i])
    #                lhs = WM._get_head_loss_oa_binary(qp[pipe_ids[i]], y, qp_vals[i], exponent)
    #                con = WM.JuMP.@build_constraint(r[i] * lhs <= inv(L[i]) * dhp[pipe_ids[i]])
    #                WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
    #            elseif abs(dhn_vals[i] - dhn_ests[i]) > 1.0e-1
    #                y = WM.var(wm, nw, :y_pipe, pipe_ids[i])
    #                lhs = WM._get_head_loss_oa_binary(qn[pipe_ids[i]], 1.0 - y, qn_vals[i], exponent)
    #                con = WM.JuMP.@build_constraint(r[i] * lhs <= inv(L[i]) * dhn[pipe_ids[i]])
    #                WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
    #            end
    #        end
    #    end

    #    for nw in WM.nw_ids(wm)
    #        pump_ids = collect(WM.ids(wm, nw, :pump))

    #        qp, g, z  = WM.var(wm, nw, :qp_pump), WM.var(wm, nw, :g_pump), WM.var(wm, nw, :z_pump)
    #        qp_vals = [max(0.0, WM.JuMP.callback_value(cb_data, qp[i])) for i in pump_ids]
    #        g_vals = [max(0.0, WM.JuMP.callback_value(cb_data, g[i])) for i in pump_ids]
    #        z_vals = [WM.JuMP.callback_value(cb_data, z[i]) >= 0.5 ? 1.0 : 0.0 for i in pump_ids]

    #        head_curves = [WM.ref(wm, nw, :pump, i)["head_curve"] for i in pump_ids]
    #        pcs = [WM._get_function_from_head_curve(head_curves[i]) for i in 1:length(pump_ids)]
    #        g_ests = [(pcs[i][1] * qp_vals[i]^2) +
    #                  (pcs[i][2] * qp_vals[i]) +
    #                  pcs[i][3] * z_vals[i] for i in 1:length(pump_ids)]

    #        for i in 1:length(pump_ids)
    #            if abs(g_vals[i] - g_ests[i]) > 1.0e-1
    #                rhs = WM._get_head_gain_oa(qp[pump_ids[i]], z[pump_ids[i]], qp_vals[i], pcs[i])
    #                con = WM.JuMP.@build_constraint(g[pump_ids[i]] <= rhs)
    #                WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
    #            end
    #        end
    #    end

    #    indicator_vars = _get_indicator_variables(wm)
    #    indicator_vals = WM.JuMP.callback_value.(Ref(cb_data), indicator_vars)
    #    is_integer = all(x -> isapprox(x, 1.0; atol = 0.01) ||
    #                     isapprox(x, 0.0; atol = 0.01), indicator_vals)

    #    if is_integer
    #        # Populate the solution of wm_nlp to use in the feasibility check.
    #        _populate_solution!(cb_data, wm, wm_nlp)
    #        wn, wnres = _simulate_solution(wm_nlp, network_path)
    #        _update_pipe_flows!(wm_nlp, wn, wnres)
    #        _update_pump_flows!(wm_nlp, wn, wnres)

    #        for nw in WM.nw_ids(wm)
    #            exponent = WM.ref(wm, nw, :alpha)
    #            nw_sol = wm_nlp.solution["nw"][string(nw)]
    #            pipe_ids = collect(WM.ids(wm, nw, :pipe))
    #            q_sol = [nw_sol["pipe"][string(i)]["q"] for i in pipe_ids]

    #            L = [WM.ref(wm, nw, :pipe, i)["length"] for i in pipe_ids]
    #            r = [WM.ref(wm, nw, :resistance, i)[1] for i in pipe_ids]
    #            qp, qn = WM.var(wm, nw, :qp_pipe), WM.var(wm, nw, :qn_pipe)
    #            dhp, dhn = WM.var(wm, nw, :dhp_pipe), WM.var(wm, nw, :dhn_pipe)

    #            for i in 1:length(pipe_ids)
    #                if q_sol[i] > 0.0
    #                    y = WM.var(wm, nw, :y_pipe, pipe_ids[i])
    #                    lhs = WM._get_head_loss_oa_binary(qp[pipe_ids[i]], y, q_sol[i], exponent)
    #                    con = WM.JuMP.@build_constraint(r[i] * lhs <= inv(L[i]) * dhp[pipe_ids[i]])
    #                    WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
    #                elseif q_sol[i] < 0.0
    #                    y = WM.var(wm, nw, :y_pipe, pipe_ids[i])
    #                    lhs = WM._get_head_loss_oa_binary(qn[pipe_ids[i]], 1.0 - y, -q_sol[i], exponent)
    #                    con = WM.JuMP.@build_constraint(r[i] * lhs <= inv(L[i]) * dhn[pipe_ids[i]])
    #                    WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
    #                end
    #            end
    #        end

    #        for nw in WM.nw_ids(wm)
    #            pump_ids = collect(WM.ids(wm, nw, :pump))
    #            qp, g, z  = WM.var(wm, nw, :qp_pump), WM.var(wm, nw, :g_pump), WM.var(wm, nw, :z_pump)
    #            head_curves = [WM.ref(wm, nw, :pump, i)["head_curve"] for i in pump_ids]
    #            pcs = [WM._get_function_from_head_curve(head_curves[i]) for i in 1:length(pump_ids)]

    #            nw_sol = wm_nlp.solution["nw"][string(nw)]
    #            q_sol = [nw_sol["pump"][string(i)]["q"] for i in pump_ids]

    #            for i in 1:length(pump_ids)
    #                if q_sol[i] > 1.0e-6
    #                    rhs = WM._get_head_gain_oa(qp[pump_ids[i]], z[pump_ids[i]], q_sol[i], pcs[i])
    #                    con = WM.JuMP.@build_constraint(g[pump_ids[i]] <= rhs)
    #                    WM._MOI.submit(wm.model, WM._MOI.UserCut(cb_data), con)
    #                end
    #            end
    #        end
    #    end
    #end

    # Register the user cut callback with the JuMP modeling object.
    #WM._MOI.set(wm.model, WM._MOI.UserCutCallback(), user_cut_callback)

    # Add additional specialized cuts.
    add_pump_volume_cuts!(wm)

    # Solve the OWF optimization problem.
    result = WM.optimize_model!(wm)

    # Save relevant algorithm metadata within the result object.
    result["objective_comparison"] = objective_comparison_table
    result["num_infeasible_solutions"] = num_infeasible_solutions

    # Return the optimization result dictionary.
    return result
end


function _populate_solution!(cb_data, wm_cb::AbstractWaterModel, wm::AbstractWaterModel)
    for nw in WM.nw_ids(wm_cb)
        for comp_type in [:pump, :regulator, :valve]
            for (i, comp) in WM.ref(wm_cb, nw, comp_type)
                var_sym = Symbol("z_" * string(comp_type))
                var = WM.var(wm_cb, nw, var_sym, i)
                val = round(WM.JuMP.callback_value(cb_data, var))
                wm.solution["nw"][string(nw)][string(comp_type)][string(i)]["status"] = val
            end
        end
    end

    nw_1 = sort(collect(WM.nw_ids(wm)))[1]

    for (i, tank) in WM.ref(wm_cb, nw_1, :tank)
        h_var = WM.var(wm_cb, nw_1, :h, tank["node"])
        h_val = WM.JuMP.callback_value(cb_data, h_var)
        p_val = h_val - WM.ref(wm_cb, nw_1, :node, tank["node"])["elevation"]
        wm.solution["nw"][string(nw_1)]["node"][string(tank["node"])]["p"] = p_val
    end
end


function _get_indicator_variables_to_nw(wm::AbstractWaterModel, nw_last::Int)
    vars = Array{WM.JuMP.VariableRef, 1}()
    network_ids = sort(collect(WM.nw_ids(wm)))[1:nw_last]

    for var_sym in [:z_pump, :z_regulator, :z_valve]
        for nw_id in network_ids
            append!(vars, vcat(WM.var(wm, nw_id, var_sym)...))
        end
    end

    return vars
end


function _get_pump_indicators(wm::AbstractWaterModel)
    vars = Array{WM.JuMP.VariableRef, 1}()

    for nw_id in WM.nw_ids(wm)
        append!(vars, vcat(WM.var(wm, nw_id, :z_pump)...))
    end

    return vars
end


function _get_indicator_variables(wm::AbstractWaterModel)
    vars = Array{WM.JuMP.VariableRef, 1}()

    for var_sym in [:z_pump, :z_regulator, :z_valve]
        for nw_id in WM.nw_ids(wm)
            append!(vars, vcat(WM.var(wm, nw_id, var_sym)...))
        end
    end

    return vars
end


function _simulate_solution(wm::AbstractWaterModel, network_path::String)
    wm_solution = Dict{String, Any}("solution" => wm.solution)
    return WMA.simulate(wm.data, wm_solution, network_path)
end


function _update_pipe_flows!(wm::AbstractWaterModel, wn, wnres)
    wm_solution = Dict{String, Any}("solution" => wm.solution)

    for i in WM.ids(wm, :pipe) # Get the list of pipe indices.
        df = WMA.get_pipe_dataframe(wm.data, wm_solution, wn, wnres, string(i))

        for nw in WM.nw_ids(wm)
            wm.solution["nw"][string(nw)]["pipe"][string(i)]["q"] = df[nw, :].flow_wntr
        end
    end
end


function _update_pump_flows!(wm::AbstractWaterModel, wn, wnres)
    wm_solution = Dict{String, Any}("solution" => wm.solution)

    for i in WM.ids(wm, :pump) # Get the list of pump indices.
        df = WMA.get_pump_dataframe(wm.data, wm_solution, wn, wnres, string(i))

        for nw in WM.nw_ids(wm)
            wm.solution["nw"][string(nw)]["pump"][string(i)]["q"] = df[nw, :].flow_wntr
        end
    end
end


function _set_van_zyl_statuses(wm::AbstractWaterModel)
    statuses = Dict{String, Array}("pmp1" => [0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1],
                                   "pmp2" => [0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 0 1 1 1 1 1 1 0],
                                   "pmp6" => [0 0 0 0 0 0 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1])

    for (nw, nw_ref) in WM.nws(wm)
        for (a, pump) in WM.ref(wm, nw, :pump)
            pump_name = pump["name"]
            pump_status = statuses[pump_name][nw]
            pump["z_min"] = Float64(pump_status)
            pump["z_max"] = Float64(pump_status)

            z_pump = WM.var(wm, nw, :z_pump, a)
            WM._fix_indicator_variable(z_pump, pump, "z")
        end
    end
end


function _calc_node_infeasibilities(wm::AbstractWaterModel, wntr_result::WMA.PyCall.PyObject)
    network_ids = sort(collect(WM.nw_ids(wm)))
    infeasibilities = Dict(i => zeros(length(network_ids)) for i in WM.ids(wm, :node))

    for (i, node) in WM.ref(wm, :node)
        node_id = string(i) # Index of the node corresponding to the node.
        head = WMA.PyCall.getproperty(wntr_result.node["head"], node_id).values[1:end-1]
        h_min = [get(WM.ref(wm, nw, :node, i), "h_min", -Inf) for nw in network_ids]
        h_max = [get(WM.ref(wm, nw, :node, i), "h_max", Inf) for nw in network_ids]
        infeasibilities[i] = min.(head .- h_min, 0.0) + max.(head .- h_max, 0.0)
    end

    return infeasibilities
end


function calc_first_nw_infeasibility(wm::AbstractWaterModel, wntr_result::WMA.PyCall.PyObject)
    node_infeasibilities = _calc_node_infeasibilities(wm, wntr_result)
    node_nw = _get_first_nw_infeasibility(wm, node_infeasibilities)

    pipe_infeasibilities = _calc_edge_infeasibilities(wm, wntr_result, :pipe)
    pipe_nw = _get_first_nw_infeasibility(wm, pipe_infeasibilities)
    pump_infeasibilities = _calc_edge_infeasibilities(wm, wntr_result, :pump)
    pump_nw = _get_first_nw_infeasibility(wm, pump_infeasibilities)

    valve_infeasibilities = _calc_edge_infeasibilities(wm, wntr_result, :valve)
    valve_nw = _get_first_nw_infeasibility(wm, valve_infeasibilities)
    short_pipe_infeasibilities = _calc_edge_infeasibilities(wm, wntr_result, :short_pipe)
    short_pipe_nw = _get_first_nw_infeasibility(wm, short_pipe_infeasibilities)
    regulator_infeasibilities = _calc_edge_infeasibilities(wm, wntr_result, :regulator)
    regulator_nw = _get_first_nw_infeasibility(wm, regulator_infeasibilities)

    min_nws = [node_nw, pipe_nw, pump_nw, valve_nw, short_pipe_nw, regulator_nw]
    filtered_array = Int64.(filter(x -> x != nothing, min_nws))
    return length(filtered_array) > 0 ? minimum(filtered_array) : nothing
end


function _get_first_nw_infeasibility(wm::AbstractWaterModel, infeasibilities::Dict{Int, Array{Float64, 1}})
    network_ids = sort(collect(WM.nw_ids(wm)))
    nw = nothing # Index of the first infeasibility.

    for (index, array) in infeasibilities
        nw_index = findfirst(x -> x > 0.0, abs.(array))
        first_nw = nw_index == nothing ? nothing : network_ids[nw_index]
        nw = first_nw == nothing ? nw : (nw == nothing ? first_nw : min(nw, first_nw))
    end

    return nw
end


function _calc_edge_infeasibilities(wm::AbstractWaterModel, wntr_result::WMA.PyCall.PyObject, sym::Symbol)
    network_ids = sort(collect(WM.nw_ids(wm)))
    infeasibilities = Dict(i => zeros(length(network_ids)) for i in WM.ids(wm, sym))

    for (i, comp) in WM.ref(wm, sym)
        comp_id = string(sym) * string(i) # String index of the pipe corresponding.
        flow = WMA.PyCall.getproperty(wntr_result.link["flowrate"], comp_id).values[1:end-1]
        q_min = [get(WM.ref(wm, nw, sym, i), "q_min", -Inf) for nw in network_ids]
        q_max = [get(WM.ref(wm, nw, sym, i), "q_max", Inf) for nw in network_ids]
        infeasibilities[i] = min.(flow .- q_min, 0.0) + max.(flow .- q_max, 0.0)
        infeasibilities[i][abs.(infeasibilities[i]) .< 1.0e-6] .= 0.0
    end

    return infeasibilities
end


function _calc_tank_infeasibilities(wm::AbstractWaterModel, wntr_result::WMA.PyCall.PyObject)
    network_ids = sort(collect(WM.nw_ids(wm)))
    infeasibilities = Dict(i => zeros(length(network_ids)) for i in WM.ids(wm, :tank))

    for (tank_id, tank) in WM.ref(wm, :tank)
        node_id = string(tank["node"]) # Index of the node corresponding to the tank.
        head = WMA.PyCall.getproperty(wntr_result.node["head"], node_id).values[1:end-1]
        elevation = [WM.ref(wm, nw, :node, tank["node"])["elevation"] for nw in network_ids]

        h_min = [get(WM.ref(wm, nw, :node, tank["node"]), "h_min", -Inf) for nw in network_ids]
        h_max = [get(WM.ref(wm, nw, :node, tank["node"]), "h_max", Inf) for nw in network_ids]

        infeasibilities[tank_id] = min.(head .- h_min, 0.0) + max.(head .- h_max, 0.0)
        infeasibilities[tank_id][end] = min(0.0, head[end] - head[1])
    end

    return infeasibilities
end


function _calc_wntr_objective(wm::AbstractWaterModel, wntr_result::WMA.PyCall.PyObject)
    wntr_objective_value, network_ids = 0.0, sort(collect(WM.nw_ids(wm)))

    for (pump_id, pump) in WM.ref(wm, :pump)
        i = "pump" * string(pump_id) # WMA pump index (a string).
        price = [WM.ref(wm, nw, :pump, pump_id)["energy_price"] for nw in network_ids]
        flow = WMA.PyCall.getproperty(wntr_result.link["flowrate"], i).values[1:end-1]
        status = WMA.PyCall.getproperty(wntr_result.link["status"], i).values[1:end-1]
        energy = WM._calc_pump_energy.(Ref(wm), Ref(wm.cnw), Ref(pump_id), Float64.(flow))
        wntr_objective_value += sum(round.(status) .* energy .* price)
    end

    return wntr_objective_value
end

function solve_owf(network_path::String, obbt_optimizer, owf_optimizer, nlp_optimizer)
    # Read in the original network data.
    network = WM.parse_file(network_path)

    # Tighten the bounds in the network.
    ext = Dict{Symbol, Any}(:pipe_breakpoints => 10, :pump_breakpoints => 10)
    WM.run_obbt_owf!(network, obbt_optimizer; model_type = LRDWaterModel, solve_relaxed = false, ext=ext)

    # Get pairwise cutting planes from the network-relaxed problem.
    wm = instantiate_model(network, LRDWaterModel, build_owf; ext=ext)
    WM.JuMP.set_optimizer(wm.model, obbt_optimizer)
    problem_sets = WM._get_pairwise_problem_sets(wm; nw = wm.cnw)
    cuts = WM._compute_pairwise_cuts!(wm, problem_sets)

    # Construct the OWF model.
    network_mn = WM.make_multinetwork(network)
    ext = Dict(:pipe_breakpoints => 10, :pump_breakpoints => 5)
    wm = WM.instantiate_model(network_mn, LRDWaterModel, build_mn_owf; ext=ext)

    # Introduce an auxiliary variable for the objective and constrain it.
    objective_function = WM.JuMP.objective_function(wm.model)
    objective_var = WM.JuMP.@variable(wm.model, base_name = "objective_auxiliary")
    WM.JuMP.@objective(wm.model, WM._MOI.MIN_SENSE, objective_var)
    WM.JuMP.@constraint(wm.model, objective_function <= objective_var)

    # Set the optimizer and other important solver parameters.
    WM.JuMP.set_optimizer(wm.model, owf_optimizer)
    WM._MOI.set(wm.model, WM._MOI.NumberOfThreads(), 1)

    # Add the pairwise cuts obtained from the relaxed problem to the OWF problem.
    for nw_id in WM.nw_ids(wm)
        # Use the same cuts for all subnetworks of the multinetwork.
        map(x -> x.variable_index_1.network_index = nw_id, cuts)
        map(x -> x.variable_index_2.network_index = nw_id, cuts)

        # Add the collection of pairwise cuts for the subnetwork.
        WM._add_pairwise_cuts!(wm, cuts)
    end

    # Solve the convex, continuously-relaxed optimal water flow problem.
    network_mn_nlp, nlp_result = deepcopy(network_mn), Dict{String, Any}()
    wm_nlp = WM.instantiate_model(network_mn_nlp, CRDWaterModel, build_mn_owf)
    WM.relax_all_binary_variables!(wm_nlp)
    nlp_result = WM.optimize_model!(wm_nlp; optimizer = nlp_optimizer)

    function lazy_cut_callback(cb_data) # Define the lazy cut callback function.
        # Populate the solution of wm_nlp to use in the feasibility check.
        _populate_solution!(cb_data, wm, wm_nlp)
        solution_is_feasible = _check_feasibility_wntr(wm_nlp, network_path)

        vars = _get_indicator_variables(wm)
        zero_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 0.0, vars)
        one_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 1.0, vars)

        if !solution_is_feasible
            # If the solution is not feasible (according to a comparison with WNTR), add a no-good cut.
            con = WM.JuMP.@build_constraint(sum(zero_vars) - sum(one_vars) >= 1.0 - length(one_vars))
            WM._MOI.submit(wm.model, WM._MOI.LazyConstraint(cb_data), con)
            true_objective = _calc_wntr_objective(wm_nlp, network_path)
        elseif solution_is_feasible
            true_objective = _calc_wntr_objective(wm_nlp, network_path)
            objective_value = WM.JuMP.callback_value(cb_data, objective_var)
            objective_difference = true_objective - objective_value

            # TODO: Add a cut based on the above below.
            con = WM.JuMP.@build_constraint(objective_var >= true_objective - true_objective * (length(one_vars) - sum(one_vars) + sum(zero_vars)))
            WM._MOI.submit(wm.model, WM._MOI.LazyConstraint(cb_data), con)
        end
    end

    # Register the lazy cut callback with the JuMP modeling object.
    WM._MOI.set(wm.model, WM._MOI.LazyConstraintCallback(), lazy_cut_callback)

    # Solve the OWF optimization problem.
    result = WM.optimize_model!(wm)
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


function _add_owf_feasibility_cut!(wm::AbstractWaterModel)
    vars = _get_indicator_variables(wm)
    zero_vars = filter(x -> round(WM.JuMP.value(x)) == 0.0, vars)
    one_vars = filter(x -> round(WM.JuMP.value(x)) == 1.0, vars)
    WM.JuMP.@constraint(wm.model, sum(zero_vars) - sum(one_vars) >= 1 - length(one_vars))
end


function _check_feasibility_wntr(wm::AbstractWaterModel, network_path::String)
    # Simulate in WNTR with the WaterModels control solution.
    wm_solution = Dict{String, Any}("solution" => wm.solution)
    wn, wnres = WMA.simulate(wm.data, wm_solution, network_path)
    is_feasible = true # This flag indicates if the control solution is feasible.

    # Check feasibility of nodal head values.
    for node_id in WM.ids(wm, :node)
        df = WMA.get_node_dataframe(wm.data, wm_solution, wn, wnres, string(node_id))

        for nw in WM.nw_ids(wm)
            node = WM.ref(wm, nw, :node, node_id)
            h_min_satisfied = haskey(node, "h_min") ? df[nw, :].head_wntr >= node["h_min"] : true
            h_max_satisfied = haskey(node, "h_max") ? df[nw, :].head_wntr <= node["h_max"] : true

            if !h_min_satisfied || !h_max_satisfied
                is_feasible = false
                break
            end
        end

        if !is_feasible
            break
        end
    end

    if !is_feasible
        return is_feasible
    end
    
    for pipe_id in WM.ids(wm, :pipe)
        df = WMA.get_pipe_dataframe(wm.data, wm_solution, wn, wnres, string(pipe_id))

        for nw in WM.nw_ids(wm)
            pipe = WM.ref(wm, nw, :pipe, pipe_id)

            if df[nw, :].flow_wntr < pipe["q_min"] || df[nw, :].flow_wntr > pipe["q_max"]
                is_feasible = false
                break
            end
        end

        if !is_feasible
            break
        end
    end

    if !is_feasible
        return is_feasible
    end

    for pump_id in WM.ids(wm, :pump)
        df = WMA.get_pump_dataframe(wm.data, wm_solution, wn, wnres, string(pump_id))

        for nw in WM.nw_ids(wm)
            pump = WM.ref(wm, nw, :pump, pump_id)

            if df[nw, :].flow_wntr > pump["q_max"]
                is_feasible = false
                break
            end
        end

        if !is_feasible
            break
        end
    end

    if !is_feasible
        return is_feasible
    end

    for short_pipe_id in WM.ids(wm, :short_pipe)
        df = WMA.get_short_pipe_dataframe(wm.data, wm_solution, wn, wnres, string(short_pipe_id))

        for nw in WM.nw_ids(wm)
            short_pipe = WM.ref(wm, nw, :short_pipe, short_pipe_id)

            if df[nw, :].flow_wntr < short_pipe["q_min"] || df[nw, :].flow_wntr > short_pipe["q_max"]
                is_feasible = false
                break
            end
        end

        if !is_feasible
            break
        end
    end

    if !is_feasible
        return is_feasible
    end

    for valve_id in WM.ids(wm, :valve)
        df = WMA.get_valve_dataframe(wm.data, wm_solution, wn, wnres, string(valve_id))

        for nw in WM.nw_ids(wm)
            valve = WM.ref(wm, nw, :valve, valve_id)

            if df[nw, :].flow_wntr < valve["q_min"] || df[nw, :].flow_wntr > valve["q_max"]
                is_feasible = false
                break
            end
        end

        if !is_feasible
            break
        end
    end

    return is_feasible
end


function _calc_wntr_objective(wm::AbstractWaterModel, network_path::String)
    wm_solution = Dict{String, Any}("solution" => wm.solution)
    wn, wnres = WMA.simulate(wm.data, wm_solution, network_path)
    wntr_objective_value = 0.0

    for pump_id in WM.ids(wm, :pump)
        df = WMA.get_pump_dataframe(wm.data, wm_solution, wn, wnres, string(pump_id))
        wntr_objective_value += sum(df.cost_wntr)
    end

    return wntr_objective_value
end

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
    ext = Dict(:pipe_breakpoints => 10, :pump_breakpoints => 10)
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
        wn, wnres = _simulate_solution(wm_nlp, network_path)

        # Calculate infeasibilities for components.
        node_infeas = _calc_node_infeasibilities(wm_nlp, wn, wnres)
        pipe_infeas = _calc_pipe_infeasibilities(wm_nlp, wn, wnres)
        pump_infeas = _calc_pump_infeasibilities(wm_nlp, wn, wnres)
        regulator_infeas = _calc_regulator_infeasibilities(wm_nlp, wn, wnres)
        short_pipe_infeas = _calc_short_pipe_infeasibilities(wm_nlp, wn, wnres)
        valve_infeas = _calc_valve_infeasibilities(wm_nlp, wn, wnres)

        # Collect the current integer solution into "zero" and "one" buckets.
        vars = _get_indicator_variables(wm) # All relevant component status variables.
        zero_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 0.0, vars)
        one_vars = filter(x -> round(WM.JuMP.callback_value(cb_data, x)) == 1.0, vars)

        if sum(node_infeas) + sum(pipe_infeas) + sum(pump_infeas) + sum(regulator_infeas) +
            sum(short_pipe_infeas) + sum(valve_infeas) > 0.0 # If any infeasibility exists...
            # If the solution is not feasible (according to a comparison with WNTR), add a no-good cut.
            con = WM.JuMP.@build_constraint(sum(zero_vars) - sum(one_vars) >= 1.0 - length(one_vars))
            WM._MOI.submit(wm.model, WM._MOI.LazyConstraint(cb_data), con)
        else
            true_objective = _calc_wntr_objective(wm_nlp, wn, wnres)
            bin_expr = true_objective * (length(one_vars) - sum(one_vars) + sum(zero_vars))
            con = WM.JuMP.@build_constraint(objective_var >= true_objective - bin_expr)
            WM._MOI.submit(wm.model, WM._MOI.LazyConstraint(cb_data), con)
        end
    end

    # Register the lazy cut callback with the JuMP modeling object.
    WM._MOI.set(wm.model, WM._MOI.LazyConstraintCallback(), lazy_cut_callback)

    # Solve the OWF optimization problem.
    return WM.optimize_model!(wm)
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
    WM.JuMP.@constraint(wm.model, sum(zero_vars) - sum(one_vars) >= 1.0 - length(one_vars))
end


function _simulate_solution(wm::AbstractWaterModel, network_path::String)
    wm_solution = Dict{String, Any}("solution" => wm.solution)
    return WMA.simulate(wm.data, wm_solution, network_path)
end


function _calc_node_infeasibilities(wm::AbstractWaterModel, wn, wnres)
    node_ids = WM.ids(wm, :node) # Get the list of node indices.
    wm_solution = Dict{String, Any}("solution" => wm.solution)
    dfs = WMA.get_node_dataframe.(Ref(wm.data), Ref(wm_solution), Ref(wn), Ref(wnres), string.(node_ids))
    inf = ones(Float64, length(WM.nw_ids(wm)), length(node_ids))

    for nw in WM.nw_ids(wm) # Loop over all subnetworks of the multinetwork.
        lbs = [get(WM.ref(wm, nw, :node, i), "h_min", -Inf) for i in node_ids]
        ubs = [get(WM.ref(wm, nw, :node, i), "h_max", Inf) for i in node_ids]
        vals = [dfs[k][nw, :].head_wntr for k in 1:length(node_ids)]
        inf[nw, :] = max.(max.(lbs .- vals, 0.0), max.(vals .- ubs, 0.0))
    end

    inf[inf .< 1.0e-3] .= 0.0 # Replace small infeasibilities with zero.
    return inf # Return the matrix of node head infeasibilities.
end


function _calc_pipe_infeasibilities(wm::AbstractWaterModel, wn, wnres)
    pipe_ids = WM.ids(wm, :pipe) # Get the list of pipe indices.
    wm_solution = Dict{String, Any}("solution" => wm.solution)
    dfs = WMA.get_pipe_dataframe.(Ref(wm.data), Ref(wm_solution), Ref(wn), Ref(wnres), string.(pipe_ids))
    inf = ones(Float64, length(WM.nw_ids(wm)), length(pipe_ids))

    for nw in WM.nw_ids(wm) # Loop over all subnetworks of the multinetwork.
        lbs = [get(WM.ref(wm, nw, :pipe, i), "q_min", -Inf) for i in pipe_ids]
        ubs = [get(WM.ref(wm, nw, :pipe, i), "q_max", Inf) for i in pipe_ids]
        vals = [dfs[k][nw, :].flow_wntr for k in 1:length(pipe_ids)]
        inf[nw, :] = max.(max.(lbs .- vals, 0.0), max.(vals .- ubs, 0.0))
    end

    inf[inf .< 1.0e-6] .= 0.0 # Replace small infeasibilities with zero.
    return inf # Return the matrix of pipe flow infeasibilities.
end


function _calc_pump_infeasibilities(wm::AbstractWaterModel, wn, wnres)
    pump_ids = WM.ids(wm, :pump) # Get the list of pump indices.
    wm_solution = Dict{String, Any}("solution" => wm.solution)
    dfs = WMA.get_pump_dataframe.(Ref(wm.data), Ref(wm_solution), Ref(wn), Ref(wnres), string.(pump_ids))
    inf = ones(Float64, length(WM.nw_ids(wm)), length(pump_ids))

    for nw in WM.nw_ids(wm) # Loop over all subnetworks of the multinetwork.
        lbs = [get(WM.ref(wm, nw, :pump, i), "q_min", -Inf) for i in pump_ids]
        ubs = [get(WM.ref(wm, nw, :pump, i), "q_max", Inf) for i in pump_ids]
        vals = [dfs[k][nw, :].flow_wntr for k in 1:length(pump_ids)]
        inf[nw, :] = max.(max.(lbs .- vals, 0.0), max.(vals .- ubs, 0.0))
    end

    inf[inf .< 1.0e-6] .= 0.0 # Replace small infeasibilities with zero.
    return inf # Return the matrix of pump flow infeasibilities.
end


function _calc_regulator_infeasibilities(wm::AbstractWaterModel, wn, wnres)
    regulator_ids = WM.ids(wm, :regulator) # Get the list of regulator indices.
    wm_solution = Dict{String, Any}("solution" => wm.solution)
    dfs = WMA.get_regulator_dataframe.(Ref(wm.data), Ref(wm_solution), Ref(wn), Ref(wnres), string.(regulator_ids))
    inf = ones(Float64, length(WM.nw_ids(wm)), length(regulator_ids))

    for nw in WM.nw_ids(wm) # Loop over all subnetworks of the multinetwork.
        lbs = [get(WM.ref(wm, nw, :regulator, i), "q_min", -Inf) for i in regulator_ids]
        ubs = [get(WM.ref(wm, nw, :regulator, i), "q_max", Inf) for i in regulator_ids]
        vals = [dfs[k][nw, :].flow_wntr for k in 1:length(regulator_ids)]
        inf[nw, :] = max.(max.(lbs .- vals, 0.0), max.(vals .- ubs, 0.0))
    end

    inf[inf .< 1.0e-6] .= 0.0 # Replace small infeasibilities with zero.
    return inf # Return the matrix of regulator flow infeasibilities.
end


function _calc_short_pipe_infeasibilities(wm::AbstractWaterModel, wn, wnres)
    short_pipe_ids = WM.ids(wm, :short_pipe) # Get the list of short_pipe indices.
    wm_solution = Dict{String, Any}("solution" => wm.solution)
    dfs = WMA.get_short_pipe_dataframe.(Ref(wm.data), Ref(wm_solution), Ref(wn), Ref(wnres), string.(short_pipe_ids))
    inf = ones(Float64, length(WM.nw_ids(wm)), length(short_pipe_ids))

    for nw in WM.nw_ids(wm) # Loop over all subnetworks of the multinetwork.
        lbs = [get(WM.ref(wm, nw, :short_pipe, i), "q_min", -Inf) for i in short_pipe_ids]
        ubs = [get(WM.ref(wm, nw, :short_pipe, i), "q_max", Inf) for i in short_pipe_ids]
        vals = [dfs[k][nw, :].flow_wntr for k in 1:length(short_pipe_ids)]
        inf[nw, :] = max.(max.(lbs .- vals, 0.0), max.(vals .- ubs, 0.0))
    end

    inf[inf .< 1.0e-6] .= 0.0 # Replace small infeasibilities with zero.
    return inf # Return the matrix of short_pipe flow infeasibilities.
end


function _calc_valve_infeasibilities(wm::AbstractWaterModel, wn, wnres)
    valve_ids = WM.ids(wm, :valve) # Get the list of valve indices.
    wm_solution = Dict{String, Any}("solution" => wm.solution)
    dfs = WMA.get_valve_dataframe.(Ref(wm.data), Ref(wm_solution), Ref(wn), Ref(wnres), string.(valve_ids))
    inf = ones(Float64, length(WM.nw_ids(wm)), length(valve_ids))

    for nw in WM.nw_ids(wm) # Loop over all subnetworks of the multinetwork.
        lbs = [get(WM.ref(wm, nw, :valve, i), "q_min", -Inf) for i in valve_ids]
        ubs = [get(WM.ref(wm, nw, :valve, i), "q_max", Inf) for i in valve_ids]
        vals = [dfs[k][nw, :].flow_wntr for k in 1:length(valve_ids)]
        inf[nw, :] = max.(max.(lbs .- vals, 0.0), max.(vals .- ubs, 0.0))
    end

    inf[inf .< 1.0e-6] .= 0.0 # Replace small infeasibilities with zero.
    return inf # Return the matrix of valve flow infeasibilities.
end


function _calc_wntr_objective(wm::AbstractWaterModel, wn, wnres)
    wm_solution = Dict{String, Any}("solution" => wm.solution)
    wntr_objective_value = 0.0

    for pump_id in [string(i) for i in WM.ids(wm, :pump)]
        df = WMA.get_pump_dataframe(wm.data, wm_solution, wn, wnres, pump_id)
        wntr_objective_value += sum(df.cost_wntr)
    end

    return wntr_objective_value
end

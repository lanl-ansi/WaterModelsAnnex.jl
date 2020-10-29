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
    WM.JuMP.set_optimizer(wm.model, owf_optimizer)

    for nw_id in WM.nw_ids(wm)
        # Use the same cuts for all subnetworks of the multinetwork.
        map(x -> x.variable_index_1.network_index = nw_id, cuts)
        map(x -> x.variable_index_2.network_index = nw_id, cuts)

        # Add the collection of pairwise cuts for the subnetwork.
        WM._add_pairwise_cuts!(wm, cuts)
    end

    # Prepare the NLP model's data dictionaries.
    network_mn_nlp, nlp_result = deepcopy(network_mn), Dict{String, Any}()

    while true
        # Solve the OWF optimization problem.
        result = WM.optimize_model!(wm)

        # Update the network data with solution data.
        WM._IM.update_data!(network_mn_nlp, result["solution"])
        WM.set_start_all!(network_mn_nlp)
        WM.fix_all_indicators!(network_mn_nlp)

        # Try to recover a feasible solution to the nonconvex problem.
        wm_nlp = WM.instantiate_model(network_mn_nlp, NCWaterModel, build_mn_owf)
        WM.relax_all_binary_variables!(wm_nlp)
        nlp_result = WM.optimize_model!(wm_nlp; optimizer = nlp_optimizer)

        # Determine whether or not to terminate the algorithm.
        if !(nlp_result["primal_status"] in [FEASIBLE_POINT, NEARLY_FEASIBLE_POINT])
            # Add the feasibility cut.
            _add_owf_feasibility_cut!(wm)
        else
            break # The solution is feasible. Terminate.
        end
    end

    return nlp_result
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

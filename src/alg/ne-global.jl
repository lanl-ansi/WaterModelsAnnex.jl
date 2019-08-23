"""
Implements an algorithm that applies feasibility cuts to network designs that
are infeasible with respect to flow and/or head bounds.
"""
function ne_global(network_path::String, modifications_path::String,
                   nlp_optimizer::JuMP.OptimizerFactory,
                   mip_optimizer::JuMP.OptimizerFactory,
                   add_feasibility_cut!::Function,
                   geojson_path::String)
    # Read in the original network data.
    network = WMs.parse_file(network_path)

    # Construct the CP feasibility model.
    wf = WMs.build_generic_model(network, WMs.NCNLPWaterModel, WMs.post_wf)
    wf_result = WaterModels.run_wf(network, NCNLPWaterModel, nlp_optimizer)

    # Add the modifications to the network data and construct the MILPR model.
    modifications = WMs.parse_file(modifications_path)
    InfrastructureModels.update_data!(network, modifications)
    InfrastructureModels.update_data!(network, wf_result["solution"])

    # Prepare warm start values for the network expansion model.
    WMs.set_start_head!(network)
    WMs.set_start_directed_head_difference!(network)
    WMs.set_start_undirected_flow_rate!(network)
    WMs.set_start_directed_flow_rate!(network)
    WMs.set_start_directed_flow_rate_ne!(network)
    WMs.set_start_flow_direction!(network)
    WMs.set_start_resistance_ne!(network)

    # Build the relaxed network expansion problem.
    ne = WMs.build_generic_model(network, WMs.MILPRWaterModel, WMs.post_ne)

    max_iterations = 10
    iteration_counter = 0
    design_is_feasible = false
    all_features = []
    println("iterations,lower_bound")

    while !design_is_feasible && iteration_counter != max_iterations
        iteration_counter += 1
        ne_solution = WMs.solve_generic_model(ne, mip_optimizer)
        objective_value = ne_solution["objective_value"]

        q, h, resistance_indices = get_cnlp_solution(ne, nlp_optimizer)
        qlb, qub, hlb, hub = check_solution_bounds(ne, q, h, resistance_indices)
        design_is_feasible = all([all(values(qlb)), all(values(qub)), all(values(hlb)), all(values(hub))])

        qlb_false = filter(x -> x.second == false, qlb)
        qub_false = filter(x -> x.second == false, qub)
        hlb_false = filter(x -> x.second == false, hlb)
        hub_false = filter(x -> x.second == false, hub)

        infeasible_arcs = sort(collect(Set(vcat(collect(keys(qlb_false)), collect(keys(qub_false))))))
        infeasible_nodes = sort(collect(Set(vcat(collect(keys(hlb_false)), collect(keys(hub_false))))))

        if !design_is_feasible
            add_feasibility_cut!(ne, nlp_optimizer)
        end

        println(iteration_counter, ",", objective_value)
    end
end

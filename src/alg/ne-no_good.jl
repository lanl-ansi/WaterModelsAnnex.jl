"""
Implements a Benders-like algorithm that applies no-good cuts to network
designs that are infeasible (with respect to flow and/or head bounds).
"""
function ne_no_good(network_path::String, modifications_path::String,
                    nlp_optimizer::JuMP.OptimizerFactory,
                    mip_optimizer::JuMP.OptimizerFactory)
    # Read in the original network data.
    network = WMs.parse_file(network_path)
    alpha = network["options"]["headloss"] == "h-w" ? 1.852 : 2.0

    # Construct the CVXNLP feasibility model.
    wf = WMs.build_generic_model(network, WMs.CVXNLPWaterModel, WMs.get_post_wf(alpha))

    # Add the modifications to the network data.
    modifications = WMs.parse_file(modifications_path)
    InfrastructureModels.update_data!(network, modifications)
    ne = WMs.build_generic_model(network, WMs.MICPWaterModel, WMs.get_post_ne(alpha))

    println(ne.model)
    optimize!(ne, mip_optimizer)

    return true
end

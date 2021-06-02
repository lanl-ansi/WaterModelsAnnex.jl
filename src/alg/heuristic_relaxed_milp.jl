function build_nw_id_partition(wm::WM.AbstractWaterModel)
    nw_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]
    return [[x] for x in Random.shuffle(nw_ids)]
end


function get_control_settings_at_nw(result::Dict{String, <:Any}, nw::Int)
    sol_nw = result["solution"]["nw"][string(nw)]
    pump_ids = sort([parse(Int, x) for x in collect(keys(sol_nw["pump"]))])
    pump_vals = [sol_nw["pump"][string(i)]["status"] for i in pump_ids]
    pump_vars = WM._VariableIndex.(nw, :pump, :z_pump, pump_ids)
    return ControlSetting(nw, pump_vars, pump_vals)
end


function get_control_settings_from_result(result::Dict{String, <:Any})
    nw_ids = sort([parse(Int, x) for x in keys(result["solution"]["nw"])])[1:end-1]
    return get_control_settings_at_nw.(Ref(result), nw_ids)
end


function feasibility_pump(wm::WM.AbstractWaterModel)
    objective = JuMP.AffExpr(0.0)

    for nw in sort(collect(WM.nw_ids(wm)))[1:end-1]
        z_pump = vcat(WM.var(wm, nw, :z_pump)...)
        z_regulator = vcat(WM.var(wm, nw, :z_regulator)...)
        z_valve = vcat(WM.var(wm, nw, :z_valve)...)

        z = vcat(z_pump, z_regulator, z_valve)
        z_one = filter(x -> round(JuMP.value(x)) == 1.0, z)
        z_zero = filter(x -> round(JuMP.value(x)) == 0.0, z)

        if length(z_one) > 0
            objective += (length(z_one) - sum(z_one))
        end

        if length(z_zero) > 0
            objective += sum(z_zero)
        end

        JuMP.unfix.(z)
        JuMP.set_binary.(z)
     end

    #  JuMP.@constraint(wm.model, objective <= 3.0)
     JuMP.@objective(wm.model, WM._MOI.MIN_SENSE, objective)
     return WM.optimize_model!(wm)
end


function feasibility_pump_at_nw(wm::WM.AbstractWaterModel, nw_partition::Array{Int64, 1}, nws_fixed::Vector{Array{Int64, 1}})
    nw_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]

    for nws_variable_counter in 1:min(length(nws_fixed), 3)
        nw_ids_current = reverse(nws_fixed)[1:nws_variable_counter]
        nw_ids_before = max.(nw_ids[1], nw_partition .- [1])
        nw_ids_after = min.(nw_ids[end], nw_partition .+ [1])
        nw_ids_variable = vcat(nw_ids_current..., nw_ids_before..., nw_ids_after...)

        # Relax variables not considered in the binary variable set.
        nw_ids_relaxed = sort(collect(setdiff(Set(nw_ids), Set(nw_ids_variable))))
        vars_to_relax = vcat(WM.get_all_binary_vars_at_nw!.(Ref(wm), nw_ids_relaxed)...)
        map(x -> JuMP.unset_binary(x), vars_to_relax)

        # Ensure variables that are being decided are not fixed.
        vars_discrete = vcat(WM.get_all_binary_vars_at_nw!.(Ref(wm), nw_ids_variable)...)
        map(x -> JuMP.is_fixed(x) && JuMP.unfix(x), vars_discrete)
        map(x -> JuMP.set_binary(x), vars_discrete)
        JuMP.optimize!(wm.model) # Solve the relaxed model.
            
        if JuMP.primal_status(wm.model) === WM._MOI.FEASIBLE_POINT
            binary_vars = filter(v -> JuMP.is_binary(v), vars_discrete)
            JuMP.fix.(binary_vars, JuMP.value.(binary_vars))
            WM.set_binary_variables!(vars_to_relax)
            return true # Found a feasible solution.
        else
            WM.set_binary_variables!(vars_to_relax)
        end
    end

    return false # Did not find a feasible solution.
end


function run_feasibility_pump(network_mn::Dict{String, Any}, mip_solver)
    # Specify model options and construct the multinetwork OWF model.
    wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf)

    # Set the optimizer and other important solver parameters.
    WM.JuMP.set_optimizer(wm.model, mip_solver)
    nw_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]
    nw_partitions = build_nw_id_partition(wm)
    
    for (i, nw_partition) in enumerate(nw_partitions)
        # Ensure a feasible solution exists for nw_partitions[1:i].
        if !feasibility_pump_at_nw(wm, nw_partition, nw_partitions[1:i])
            return nothing
        end
    end

    return WM.optimize_model!(wm)
end


function run_heuristic(network::Dict{String, Any}, mip_optimizer, nlp_optimizer)
    network_mn = WM.make_multinetwork(network)
    wm_micp = construct_owf_model_relaxed(network_mn, nlp_optimizer)
    result_micp = WM.optimize_model!(wm_micp; relax_integrality = true)
    update_multinetwork_heuristic_breakpoints!(network_mn, result_micp)
    networks_mn = [deepcopy(network_mn) for i in 1:Threads.nthreads()]

    networks = [deepcopy(network) for i in 1:Threads.nthreads()]
    heuristic_result = Vector{Any}([nothing for i in 1:Threads.nthreads()])
    wms = [_instantiate_cq_model(networks[i], nlp_optimizer) for i in 1:Threads.nthreads()]
    settings = [[ControlSetting(n, [], [])] for n in 1:Threads.nthreads()]
    results = [[SimulationResult(false, Dict{Int, Float64}(), 0.0)] for n in 1:Threads.nthreads()]
    feasible, costs = [false for n in 1:Threads.nthreads()], zeros(Float64, Threads.nthreads())
    
    Threads.@threads for k in 1:Threads.nthreads()
        heuristic_result[k] = run_feasibility_pump(networks_mn[k], mip_optimizer)
        is_feasible = heuristic_result[k] === nothing ? false : heuristic_result[k]["primal_status"] === WM._MOI.FEASIBLE_POINT
        settings[k] = is_feasible ? get_control_settings_from_result(heuristic_result[k]) : settings[k]
        results[k] = is_feasible ? simulate_control_setting.(Ref(wms[k]), settings[k]) : results[k]
        feasible[k] = is_feasible ? all(x.feasible for x in results[k]) : feasible[k]
        # costs[k] = is_feasible ? sum(x.cost for x in results[k]) : Inf
        costs[k] = is_feasible ? heuristic_result[k]["objective"] : Inf
    end

    println(minimum(costs), " ", costs)
    # feasible_ids = findall(x -> x == true, feasible)
    # println(minimum(costs[feasible_ids]), " ", costs[feasible_ids])

    # return feasible, costs
end
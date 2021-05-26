function build_nw_id_partition(wm::WM.AbstractWaterModel)
    return Random.shuffle([[x] for x in sort(collect(WM.nw_ids(wm)))[1:end-1]])

    nw_ids = [x for x in sort(collect(WM.nw_ids(wm)))[1:end-1]]
    nw_ids_relax, nw_mid = [], Int(nw_ids[Int(length(nw_ids) / 2)])

    for i in 1:Int(floor(length(nw_ids) / 2))
        nw_1, nw_2 = nw_ids[i], nw_ids[end-i+1]
        array_append = nw_1 != nw_2 ? [nw_1, nw_2] : [nw_1]
        push!(nw_ids_relax, array_append)
    end

    nw_ids_final = deepcopy(nw_ids_relax)

    for i in 1:Int(floor(length(nw_ids_relax) / 2))
        nw_ids_final[2*(i-1)+1] = nw_ids_relax[i]
        nw_ids_final[2*(i-1)+2] = nw_ids_relax[end-i+1]
    end

    # return nw_ids_final

    # return [[x] for x in vcat(nw_ids_final...)]
    return Random.shuffle([[x] for x in vcat(nw_ids_final...)])
end


function set_heads_at_tanks_nws(wm::WM.AbstractWaterModel, nw_ids::Array{Int64, 1})
    for nw in nw_ids
        for (i, tank) in WM.ref(wm, nw, :tank)
            h_var = WM.var(wm, nw, :h, tank["node"])
            h_val = JuMP.value(h_var)
            h_lb, h_ub = JuMP.lower_bound(h_var), JuMP.upper_bound(h_var)

            JuMP.set_lower_bound(h_var, max(h_lb, h_val))
            JuMP.set_upper_bound(h_var, min(h_ub, h_val))
            JuMP.fix(h_var, JuMP.value(h_var); force = true)
        end
    end
end


function create_control_settings_from_result(result::Dict{String, <:Any})
    control_settings = Array{ControlSetting, 1}([])

    for nw in sort([parse(Int, x) for x in keys(result["solution"]["nw"])])[1:end-1]
        sol_nw = result["solution"]["nw"][string(nw)]
        pump_ids = sort([parse(Int, x) for x in collect(keys(sol_nw["pump"]))])
        pump_vals = [sol_nw["pump"][string(i)]["status"] for i in pump_ids]
        pump_vars = WM._VariableIndex.(nw, :pump, :z_pump, pump_ids)
        push!(control_settings, ControlSetting(nw, pump_vars, pump_vals))
    end

    return control_settings
end


function compute_heuristic_schedule(network::Dict{String, Any}, mip_solver)
    network_mn = WM.make_multinetwork(network)

    # Specify model options and construct the multinetwork OWF model.
    ext = Dict(:pipe_breakpoints => 5, :pump_breakpoints => 5)
    wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf; ext = ext)

    # Set the optimizer and other important solver parameters.
    WM.JuMP.set_optimizer(wm.model, mip_solver)
    nw_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]
    partitions = build_nw_id_partition(wm)

    num_partitions = Int(floor(0.5 * length(partitions)))

    for nw_partition in build_nw_id_partition(wm)[1:num_partitions]
        nws_to_relax = sort(collect(setdiff(Set(nw_ids), Set(nw_partition))))
        vars_unset = vcat(WM.relax_all_binary_variables_at_nw!.(Ref(wm), nws_to_relax)...)
        JuMP.optimize!(wm.model) # Solve the relaxed model.

        if JuMP.primal_status(wm.model) !== WM._MOI.FEASIBLE_POINT
            return nothing
        end

        binary_vars = filter(v -> JuMP.is_binary(v), JuMP.all_variables(wm.model))
        binary_vars_ny = filter(v -> !occursin("_x", JuMP.name(v)), binary_vars)
        # binary_vars_nxy = filter(v -> !occursin("_y", JuMP.name(v)), binary_vars_ny)
        
        # set_heads_at_tanks_nws(wm, nw_partition .+ 1)
        JuMP.fix.(binary_vars_ny, JuMP.value.(binary_vars_ny))
        WM.set_binary_variables!(vars_unset)
    end

    result = WM.optimize_model!(wm)

    if JuMP.primal_status(wm.model) === WM._MOI.FEASIBLE_POINT
        return create_control_settings_from_result(result)
    else
        return nothing
    end
end


function run_heuristic(network::Dict{String, Any}, mip_optimizer, nlp_optimizer)
    networks = [deepcopy(network) for i in 1:Threads.nthreads()]
    wms = [_instantiate_cq_model(networks[i], nlp_optimizer) for i in 1:Threads.nthreads()]
    settings = [[ControlSetting(n, [], [])] for n in 1:Threads.nthreads()]
    results = [[SimulationResult(false, Dict{Int, Float64}(), 0.0)] for n in 1:Threads.nthreads()]
    feasible, costs = [false for n in 1:Threads.nthreads()], zeros(Float64, Threads.nthreads())
    
    Threads.@threads for k in 1:Threads.nthreads()
        heuristic_result = compute_heuristic_schedule(networks[k], mip_optimizer)
        settings[k] = heuristic_result !== nothing ? heuristic_result : settings[k]
        results[k] = heuristic_result !== nothing ? simulate_control_setting.(Ref(wms[k]), settings[k]) : results[k]
        feasible[k] = heuristic_result !== nothing ? all(x.feasible for x in results[k]) : feasible[k]
        costs[k] = heuristic_result !== nothing ? sum(x.cost for x in results[k]) : Inf
    end

    feasible_ids = findall(x -> x == true, feasible)
    println(minimum(costs[feasible_ids]), " ", costs[feasible_ids])

    return feasible, costs
end
function build_nw_id_partition(wm::WM.AbstractWaterModel; use_partition::Bool = false)
    nw_ids = [x for x in sort(collect(WM.nw_ids(wm)))[1:end-1]]
    return Random.shuffle([[x] for x in nw_ids])
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


function compute_heuristic_schedule(network::Dict{String, Any}, mip_solver; use_partition::Bool = false)
    network_mn = WM.make_multinetwork(network)

    # Specify model options and construct the multinetwork OWF model.
    ext = Dict(:pipe_breakpoints => 5, :pump_breakpoints => 5)
    wm = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf; ext = ext)

    # Set the optimizer and other important solver parameters.
    WM.JuMP.set_optimizer(wm.model, mip_solver)
    nw_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]
    
    for nw_partition in build_nw_id_partition(wm)
        nws_before = max.(nw_ids[1], nw_partition .- [1, 2, 3])
        nws_after = min.(nw_ids[end], nw_partition .+ [1, 2, 3])

        vars_before = vcat(WM.get_all_binary_vars_at_nw!.(Ref(wm), nws_before)...)
        filter(v -> occursin("_x", JuMP.name(v)) || occursin("_y", JuMP.name(v)), vars_before)

        vars_after = vcat(WM.get_all_binary_vars_at_nw!.(Ref(wm), nws_after)...)
        filter(v -> occursin("_x", JuMP.name(v)) || occursin("_y", JuMP.name(v)), vars_after)

        vars_current = vcat(WM.get_all_binary_vars_at_nw!.(Ref(wm), nw_partition)...)
        nws_to_relax = sort(collect(setdiff(Set(nw_ids), Set(nw_partition))))

        vars_to_relax = vcat(WM.get_all_binary_vars_at_nw!.(Ref(wm), nws_to_relax)...)
        # vars_to_relax = filter(v -> !occursin("_y", JuMP.name(v)), vars_to_relax)
        # vars_to_relax = filter(v -> !occursin("_x", JuMP.name(v)), vars_to_relax)

        map(x -> JuMP.unset_binary(x), vars_to_relax)
        WM.set_binary_variables!(vars_before)
        WM.set_binary_variables!(vars_after)

        JuMP.optimize!(wm.model) # Solve the relaxed model.

        if JuMP.primal_status(wm.model) !== WM._MOI.FEASIBLE_POINT
            return nothing
        end

        binary_vars = filter(v -> JuMP.is_binary(v), vars_current)
        binary_vars_ny = filter(v -> !occursin("_x", JuMP.name(v)), binary_vars)
        binary_vars_nyx = filter(v -> !occursin("_y", JuMP.name(v)), binary_vars_ny)
        
        JuMP.fix.(binary_vars_nyx, JuMP.value.(binary_vars_nyx))
        WM.set_binary_variables!(vars_to_relax)
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
        heuristic_result = compute_heuristic_schedule(networks[k], mip_optimizer; use_partition = k == 1)
        settings[k] = heuristic_result !== nothing ? heuristic_result : settings[k]
        results[k] = heuristic_result !== nothing ? simulate_control_setting.(Ref(wms[k]), settings[k]) : results[k]
        feasible[k] = heuristic_result !== nothing ? all(x.feasible for x in results[k]) : feasible[k]
        costs[k] = heuristic_result !== nothing ? sum(x.cost for x in results[k]) : Inf
    end

    feasible_ids = findall(x -> x == true, feasible)
    println(minimum(costs[feasible_ids]), " ", costs[feasible_ids])

    return feasible, costs
end
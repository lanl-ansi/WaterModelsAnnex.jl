function relax_indicators_at_nw!(wm::WM.AbstractWaterModel, nw::Int)
    var_symbols = Array{Symbol}([:z_pump, :z_regulator, :z_valve])
    vars = vcat([vcat(WM.var(wm, nw, s)...) for s in var_symbols]...)
    WM._relax_binary_variable!.(vars)
end


function unrelax_indicators_at_nw!(wm::WM.AbstractWaterModel, nw::Int)
    var_symbols = Array{Symbol}([:z_pump, :z_regulator, :z_valve])
    vars = vcat([vcat(WM.var(wm, nw, s)...) for s in var_symbols]...)
end


function run_pwlrd_heuristic(network_mn::Dict{String, <:Any}, optimizer, time_limit::Float64)
    WM.set_flow_partitions_si!(network_mn, 1.0, 1.0e-4)
    wm = WM.instantiate_model(network_mn, WM.NCWaterModel, WM.build_mn_owf)
    JuMP.set_optimizer(wm.model, optimizer)
    undo_relax = JuMP.relax_integrality(wm.model)

    nw_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]
    z_sum = sum(sum(WM.var(wm, nw, :z_pump)) for nw in nw_ids)

    JuMP.@objective(wm.model, WM.JuMP.MOI.MIN_SENSE, z_sum)
    result_relaxed_min = WM.optimize_model!(wm)
    JuMP.@objective(wm.model, WM.JuMP.MOI.MAX_SENSE, z_sum)
    result_relaxed_max = WM.optimize_model!(wm)
    undo_relax() # Undo all continuous relaxations.

    for nw in nw_ids
        sol_min = result_relaxed_min["solution"]["nw"][string(nw)]
        sol_max = result_relaxed_max["solution"]["nw"][string(nw)]

        for i in WM.ids(wm, nw, :pump)
            z = WM.var(wm, nw, :z_pump, i)
            status_min = sol_min["pump"][string(i)]["status"]
            status_max = sol_max["pump"][string(i)]["status"]

            if status_min > 0.75
                JuMP.fix(z, 1.0; force = true)
            elseif status_max < 1.0e-4
                JuMP.fix(z, 0.0; force = true)
            end
        end
    end

    WM.objective_owf(wm)
    result = WM.optimize_model!(wm)

    return result_relaxed_min
end


function get_max_z_index(result_mn::Dict{String, <:Any}, nw_ids::Vector{Int}, ids_picked::Vector{Tuple})
    z_max, z_id, nw_id = 0.0, 0, 0
    
    for nw in nw_ids
        pump_sol = result_mn["solution"]["nw"][string(nw)]["pump"]
        ids = sort(collect(keys(pump_sol)))
        vals = [pump_sol[i]["status"] for i in ids]
        best_z, max_id = findmax(vals)

        if vals[max_id] > z_max && !((nw, ids[max_id]) in ids_picked)
            z_max, z_id = vals[max_id], ids[max_id]
            nw_id = nw
        end
    end

    return nw_id, z_id
end


function run_nc_heuristic(network_mn::Dict{String, <:Any}, mip_optimizer, nlp_optimizer, time_limit::Float64)
    wm_nlp = WM.instantiate_model(network_mn, CDWaterModel, WM.build_mn_owf)
    constraint_strong_duality(wm_nlp) # Minimize this.
    JuMP.set_optimizer(wm_nlp.model, nlp_optimizer)
    undo_relax = JuMP.relax_integrality(wm_nlp.model)
    found_solution, num_iterations, ids_fixed = false, 0, Vector{Tuple}([])
    nw_ids = sort(collect(WM.nw_ids(wm_nlp)))[1:end-1]

    # z_sum = sum(sum(WM.var(wm_nlp, nw, :z_pump)) for nw in nw_ids)

    while !found_solution && num_iterations < 15
        # JuMP.@objective(wm_nlp.model, WM.JuMP.MOI.MIN_SENSE, z_sum)
        result_relaxed = WM.optimize_model!(wm_nlp)
        println(JuMP.termination_status(wm_nlp.model))

        nw_id, pump_key = get_max_z_index(result_relaxed, nw_ids, ids_fixed)
        push!(ids_fixed, (nw_id, pump_key))

        pump_sol = result_relaxed["solution"]["nw"][string(nw_id)]["pump"][pump_key]

        z_var = WM.var(wm_nlp, nw_id, :z_pump, parse(Int, pump_key))
        JuMP.fix(z_var, 1.0; force = true)
        num_iterations += 1

        println(ids_fixed)
    end

    # wm_mip = WM.instantiate_model(network_mn, WM.PWLRDWaterModel, WM.build_mn_owf)
    # JuMP.set_optimizer(wm_mip.model, mip_optimizer)

    
    

    

    # max_z, min_z = 0.0, Inf

    # for nw in nw_ids
    #     sol_min = result_relaxed_min["solution"]["nw"][string(nw)]
    #     sol_max = result_relaxed_max["solution"]["nw"][string(nw)]


    # JuMP.@objective(wm_nlp.model, WM.JuMP.MOI.MAX_SENSE, z_sum)
    # result_relaxed_max = WM.optimize_model!(wm_nlp)

    # for nw in nw_ids
    #     sol_min = result_relaxed_min["solution"]["nw"][string(nw)]
    #     sol_max = result_relaxed_max["solution"]["nw"][string(nw)]

    #     for i in WM.ids(wm_nlp, nw, :pump)
    #         z_nlp = WM.var(wm_nlp, nw, :z_pump, i)
    #         z_mip = WM.var(wm_mip, nw, :z_pump, i)
    #         status_min = sol_min["pump"][string(i)]["status"]
    #         status_max = sol_max["pump"][string(i)]["status"]

    #         if status_min > 0.80
    #             JuMP.fix(z_nlp, 1.0; force = true)
    #             JuMP.fix(z_mip, 1.0; force = true)
    #         elseif status_max < 1.0e-4
    #             JuMP.fix(z_nlp, 0.0; force = true)
    #             JuMP.fix(z_mip, 0.0; force = true)
    #         end
    #     end
    # end

    # WM.objective_owf(wm_mip)
    # return WM.optimize_model!(wm_mip)
end


function run_minlp_heuristic(wm::WM.AbstractWaterModel, time_limit::Float64)
    nw_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]
    relax_indicators_at_nw!.(Ref(wm), nw_ids[1:end])
    result_relaxation = WM.optimize_model!(wm)

    # for nw_id in reverse(nw_ids)
    #     pump_sol_nw = result_relaxation["solution"]["nw"][string(nw_id)]["pump"]
    #     pump_keys = [string(i) for i in WM.ids(wm, nw_id, :pump)]
    #     pump_statuses = [pump_sol_nw[i]["status"] for i in pump_keys]
    #     pump_ids = [parse(Int, i) for i in pump_keys]
    #     pump_vars = WM.var.(Ref(wm), nw_id, :z_pump, pump_ids)

    #     JuMP.fix.(pump_vars, round.(pump_statuses); force = true)
    #     result_relaxation = WM.optimize_model!(wm)

    #     # println(pump_vars)
    # end
end
function solve_heuristic_linear_program(wm::WM.AbstractWaterModel, control_settings::Array{ControlSetting, 1}, simulation_results::Array{SimulationResult, 1}, optimizer)
    model = JuMP.Model(optimizer)
    h = heuristic_linear_program_head_variables(wm, model)
    z = heuristic_linear_program_weight_variables(model, simulation_results)
    obj = add_heuristic_linear_program_objective!(model, z, simulation_results)
    add_heuristic_linear_program_constraints!(wm, model, h, z, control_settings, simulation_results)
    return optimize_heuristic_linear_program!(model, z)
end


function heuristic_linear_program_head_variables(wm::WM.AbstractWaterModel, model::JuMP.Model)
    head = Dict{Int, Any}(nw => JuMP.@variable(model, [i in [x["node"]
        for (i, x) in WM.ref(wm, nw, :tank)]]) for nw in WM.nw_ids(wm))
        
    for nw in WM.nw_ids(wm)
        for (i, tank) in WM.ref(wm, nw, :tank)
            node = WM.ref(wm, nw, :node, tank["node"])
            JuMP.set_lower_bound(head[nw][tank["node"]], node["head_min"])
            JuMP.set_upper_bound(head[nw][tank["node"]], node["head_max"])
        end
    end

    return head
end


function heuristic_linear_program_weight_variables(model::JuMP.Model, simulation_results::Array{SimulationResult, 1})
    return JuMP.@variable(model, [i in 1:length(simulation_results)],
        lower_bound = 0.0, upper_bound = 1.0)
end


function add_heuristic_linear_program_objective!(model::JuMP.Model, z, simulation_results::Array{SimulationResult, 1})
    cost = sum(z[i] * simulation_results[i].cost for i in 1:length(simulation_results))
    return JuMP.@objective(model, WM._MOI.MIN_SENSE, cost)
end


function add_heuristic_linear_program_constraints!(wm, model, h, z, control_settings, simulation_results)
    for nw in sort(collect(WM.nw_ids(wm)))[1:end-1]
        time_step = WM.ref(wm, nw, :time_step)
        indices = findall(x -> x.network_id == nw, control_settings)
        @assert length(indices) > 0 # There must be at least one setting per time step.
        
        for (i, tank) in WM.ref(wm, nw, :tank)
            h_cur, h_next = h[nw][tank["node"]], h[nw+1][tank["node"]]
            q_tank = sum(simulation_results[k].q_tank[i] * z[k] for k in indices)
            coeff = time_step / (0.25 * pi * tank["diameter"]^2)
            JuMP.@constraint(model, h_cur - h_next == coeff * q_tank)
        end
        
        JuMP.@constraint(model, sum(z[k] for k in indices) == 1.0)
    end
end


function optimize_heuristic_linear_program!(model::JuMP.Model, z)
    JuMP.optimize!(model)
    return JuMP.value.(z)
end
"Create flow-related variables common to all directed flow models for node-connecting components."
function WM.variable_flow(wm::AbstractCQModel; nw::Int = WM._IM.nw_id_default, bounded::Bool = false, report::Bool = true)
    for name in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
        # Create directed flow (`qp` and `qn`) variables for each component.
        WM._variable_component_flow(wm, name; nw = nw, bounded = bounded, report = report)
    end

    JuMP.fix.(WM.var(wm, nw, :qn_pump), 0.0; force = true)
    JuMP.fix.(WM.var(wm, nw, :qn_regulator), 0.0; force = true)
end


function WM.constraint_flow_conservation(
    wm::AbstractCQModel, n::Int, i::Int, pipe_fr::Vector{Int64},
    pipe_to::Vector{Int64}, des_pipe_fr::Vector{Int64}, des_pipe_to::Vector{Int64},
    pump_fr::Vector{Int64}, pump_to::Vector{Int64}, regulator_fr::Vector{Int64},
    regulator_to::Vector{Int64}, short_pipe_fr::Vector{Int64},
    short_pipe_to::Vector{Int64}, valve_fr::Vector{Int64}, valve_to::Vector{Int64},
    reservoirs::Vector{Int64}, tanks::Vector{Int64}, dispatchable_demands::Vector{Int64},
    fixed_demand::Float64)
    # Collect flow variable references per component.
    qp_pipe, qn_pipe = WM.var(wm, n, :qp_pipe), WM.var(wm, n, :qn_pipe)
    qp_pump, qn_pump = WM.var(wm, n, :qp_pump), WM.var(wm, n, :qn_pump)
    qp_regulator, qn_regulator = WM.var(wm, n, :qp_regulator), WM.var(wm, n, :qn_regulator)
    qp_short_pipe, qn_short_pipe = WM.var(wm, n, :qp_short_pipe), WM.var(wm, n, :qn_short_pipe)
    qp_valve, qn_valve = WM.var(wm, n, :qp_valve), WM.var(wm, n, :qn_valve)
    q_reservoir, q_tank = WM.var(wm, n, :q_reservoir), WM.var(wm, n, :q_tank)
    q_demand = WM.var(wm, n, :q_demand)
    
    if length(pump_fr) > 0 || length(pump_to) > 0
        z_pump = length(qp_pump) > 0 ? WM.var(wm, n, :z_pump) : 0.0
        q_pump_fr = JuMP.@NLexpression(wm.model, sum(z_pump[a] * (qp_pump[a] - qn_pump[a]) for a in pump_fr))
        q_pump_to = JuMP.@NLexpression(wm.model, sum(z_pump[a] * (qp_pump[a] - qn_pump[a]) for a in pump_to))
    else
        q_pump_fr = q_pump_to = 0.0
    end

    if length(regulator_fr) > 0 || length(regulator_to) > 0
        z_regulator = length(qp_regulator) > 0 ? WM.var(wm, n, :z_regulator) : 0.0
        q_regulator_fr = JuMP.@NLexpression(wm.model, sum(z_regulator[a] * (qp_regulator[a]) for a in regulator_fr))
        q_regulator_to = JuMP.@NLexpression(wm.model, sum(z_regulator[a] * (qp_regulator[a]) for a in regulator_to))
    else
        q_regulator_fr = q_regulator_to = 0.0
    end

    if length(valve_fr) > 0 || length(valve_to) > 0
        z_valve = length(qp_valve) > 0 ? WM.var(wm, n, :z_valve) : 0.0
        q_valve_fr = JuMP.@NLexpression(wm.model, sum(z_valve[a] * (qp_valve[a] - qn_valve[a]) for a in valve_fr))
        q_valve_to = JuMP.@NLexpression(wm.model, sum(z_valve[a] * (qp_valve[a] - qn_valve[a]) for a in valve_to))
    else
        q_valve_fr = q_valve_to = 0.0
    end
   
    if !isapprox(fixed_demand, 0.0; atol = 1.0e-7)
        fixed_demand_val = wm.model[:fixed_demands][i]
    else
        fixed_demand_val = 0.0
    end

    # Add the flow conservation constraint.
    WM.con(wm, n, :flow_conservation)[i] = JuMP.@NLconstraint(wm.model, -
        sum(qp_pipe[a] - qn_pipe[a] for a in pipe_fr) +
        sum(qp_pipe[a] - qn_pipe[a] for a in pipe_to) -
        q_pump_fr + q_pump_to - q_regulator_fr + q_regulator_to - q_valve_fr + q_valve_to -
        sum(qp_short_pipe[a] - qn_short_pipe[a] for a in short_pipe_fr) +
        sum(qp_short_pipe[a] - qn_short_pipe[a] for a in short_pipe_to) <= -
        sum(q_reservoir[k] for k in reservoirs) - sum(q_tank[k] for k in tanks) +
        sum(q_demand[k] for k in dispatchable_demands) + fixed_demand_val)

    WM.con(wm, n, :flow_conservation)[i] = JuMP.@NLconstraint(wm.model, -
        sum(qp_pipe[a] - qn_pipe[a] for a in pipe_fr) +
        sum(qp_pipe[a] - qn_pipe[a] for a in pipe_to) -
        q_pump_fr + q_pump_to - q_regulator_fr + q_regulator_to - q_valve_fr + q_valve_to -
        sum(qp_short_pipe[a] - qn_short_pipe[a] for a in short_pipe_fr) +
        sum(qp_short_pipe[a] - qn_short_pipe[a] for a in short_pipe_to) >= -
        sum(q_reservoir[k] for k in reservoirs) - sum(q_tank[k] for k in tanks) +
        sum(q_demand[k] for k in dispatchable_demands) + fixed_demand_val)
end


""
function parameter_pump_indicator(wm::AbstractCQModel; nw::Int = WM.nw_id_default)
    JuMP.@NLparameter(wm.model, z_pump[a in WM.ids(wm, nw, :pump)] == 0.0)
    WM.var(wm, nw)[:z_pump] = wm.model[:z_pump] = z_pump
end


""
function update_pump_indicator(wm::AbstractCQModel, data::Dict{Int, Float64}; nw::Int = WM.nw_id_default)
    for (i, pump) in WM.ref(wm, nw, :pump)
        JuMP.set_value(wm.model[:z_pump][i], data[i])
    end
end


""
function parameter_regulator_indicator(wm::AbstractCQModel; nw::Int = WM.nw_id_default)
    JuMP.@NLparameter(wm.model, z_regulator[a in WM.ids(wm, nw, :regulator)] == 0.0)
    WM.var(wm, nw)[:z_regulator] = wm.model[:z_regulator] = z_regulator
end


""
function update_regulator_indicator(wm::AbstractCQModel, data::Dict{Int, Float64}; nw::Int = WM.nw_id_default)
    for (i, regulator) in WM.ref(wm, nw, :regulator)
        JuMP.set_value(wm.model[:z_regulator][i], data[i])
    end
end


""
function parameter_valve_indicator(wm::AbstractCQModel; nw::Int = WM.nw_id_default)
    JuMP.@NLparameter(wm.model, z_valve[a in WM.ids(wm, nw, :valve)] == 0.0)
    WM.var(wm, nw)[:z_valve] = wm.model[:z_valve] = z_valve
end


""
function update_valve_indicator(wm::AbstractCQModel, data::Dict{Int, Float64}; nw::Int = WM.nw_id_default)
    for (i, valve) in WM.ref(wm, nw, :valve)
        JuMP.set_value(wm.model[:z_valve][i], data[i])
    end
end


""
function parameter_tank_head(wm::AbstractCQModel; nw::Int = WM.nw_id_default)
    node_ids = [tank["node"] for (i, tank) in WM.ref(wm, nw, :tank)]
    JuMP.@NLparameter(wm.model, h_tank[i in node_ids] == 0.0)
    WM.var(wm, nw)[:h_tank] = wm.model[:h_tank] = h_tank
end


""
function update_tank_head(wm::AbstractCQModel, data::Dict{Int, Float64}; nw::Int = WM.nw_id_default)
    for (i, tank) in WM.ref(wm, nw, :tank)
        node = WM.ref(wm, nw, :node, tank["node"])
        head = node["elevation"] + tank["init_level"]
        JuMP.set_value(wm.model[:h_tank][tank["node"]], head)
    end
end


""
function parameter_reservoir_head(wm::AbstractCQModel; nw::Int = WM.nw_id_default)
    node_ids = [res["node"] for (i, res) in WM.ref(wm, nw, :reservoir)]
    JuMP.@NLparameter(wm.model, h_reservoir[i in node_ids] == 0.0)
    WM.var(wm, nw)[:h_reservoir] = wm.model[:h_reservoir] = h_reservoir
end


""
function update_reservoir_head(wm::AbstractCQModel, data::Dict{Int, Float64}; nw::Int = WM.nw_id_default)
    for (i, reservoir) in WM.ref(wm, nw, :reservoir)
        node = WM.ref(wm, nw, :node, reservoir["node"])
        JuMP.set_value(wm.model[:h_reservoir][node["index"]], node["head_nominal"])
    end
end


""
function parameter_fixed_demand(wm::AbstractCQModel; nw::Int = WM.nw_id_default)
    node_ids = [demand["node"] for (i, demand) in WM.ref(wm, nw, :demand)]
    JuMP.@NLparameter(wm.model, fixed_demands[i in node_ids] == 0.0)
    WM.var(wm, nw)[:fixed_demands] = wm.model[:fixed_demands] = fixed_demands
end


""
function update_fixed_demand(wm::AbstractCQModel, data::Dict{Int, Float64}; nw::Int = WM.nw_id_default)
    for (i, demand) in WM.ref(wm, nw, :demand)
        node = WM.ref(wm, nw, :node, demand["node"])
        JuMP.set_value(wm.model[:fixed_demands][node["index"]], demand["flow_nominal"])
    end
end


""
function update_parameters(wm::AbstractCQModel, data::Dict{String, <:Any}; nw::Int = WM.nw_id_default)
    update_pump_indicator(wm, data["pump"]; nw = nw)
    update_regulator_indicator(wm, data["regulator"]; nw = nw)
    update_valve_indicator(wm, data["valve"]; nw = nw)
    update_tank_head(wm, data["tank"]; nw = nw)
    update_reservoir_head(wm, data["reservoir"]; nw = nw)
end


function objective_strong_duality(wm::AbstractCQModel)
    wm_data = WM.get_wm_data(wm.data)
    base_length = get(wm_data, "base_length", 1.0)
    base_mass = get(wm_data, "base_mass", 1.0)
    base_time = get(wm_data, "base_time", 1.0)
    
    pipe_type = wm.ref[:it][WM.wm_it_sym][:head_loss]
    viscosity = wm.ref[:it][WM.wm_it_sym][:viscosity]

    # Get NLparameter references.
    h_tank = wm.model[:h_tank]
    h_reservoir = wm.model[:h_reservoir]
    z_pump = wm.model[:z_pump]

    f_1 = Array{Any, 1}([0.0])
    f_2 = Array{Any, 1}([0.0])

    for n in WM.nw_ids(wm)
        # Get pipe flow variables.
        qp_pipe = WM.var(wm, n, :qp_pipe)
        qn_pipe = WM.var(wm, n, :qn_pipe)
        
        for (a, pipe) in WM.ref(wm, n, :pipe)
            L_x_r = pipe["length"] * WM._calc_pipe_resistance(
                pipe, pipe_type, viscosity, base_length, base_mass, base_time)
            
            push!(f_1, JuMP.@NLexpression(
                wm.model, L_x_r * head_loss(qp_pipe[a])))
            push!(f_1, JuMP.@NLexpression(
                wm.model, L_x_r * head_loss(qn_pipe[a])))
        end

        # Get design pipe flow variables.
        qp_des_pipe = WM.var(wm, n, :qp_des_pipe)
        qn_des_pipe = WM.var(wm, n, :qn_des_pipe)

        for (a, des_pipe) in WM.ref(wm, n, :des_pipe)
            L_x_r = des_pipe["length"] * WM._calc_pipe_resistance(
                des_pipe, pipe_type, viscosity, base_length, base_time)
            
            push!(f_1, JuMP.@NLexpression(
                wm.model, L_x_r * head_loss(qp_des_pipe[a])))
            push!(f_1, JuMP.@NLexpression(
                wm.model, L_x_r * head_loss(qn_des_pipe[a])))
        end

        qp_pump = WM.var(wm, n, :qp_pump)
        
        for (a, pump) in WM.ref(wm, n, :pump)
            @assert pump["pump_type"] in [WM.PUMP_QUADRATIC,
                WM.PUMP_BEST_EFFICIENCY_POINT, WM.PUMP_LINEAR_POWER]

            if pump["pump_type"] in [WM.PUMP_QUADRATIC, WM.PUMP_BEST_EFFICIENCY_POINT]
                c = WM._calc_head_curve_coefficients(pump)
            
                push!(f_2, JuMP.@NLexpression(wm.model,
                    z_pump[a] * c[1] / 3.0 * qp_pump[a]^3 +
                    z_pump[a] * 0.5 * c[2] * qp_pump[a]^2 +
                    z_pump[a] * c[3] * qp_pump[a]))
            elseif pump["pump_type"] == WM.PUMP_LINEAR_POWER
                c = WM._calc_head_curve_coefficients(pump)
                
                push!(f_2, JuMP.@NLexpression(wm.model,
                    z_pump[a] * c[1] * qp_pump[a] +
                    z_pump[a] * (1.0 / (1.0 + c[3])) *
                    c[2] * qp_pump[a]^(c[3] + 1.0)))
            end
        end

        q_tank = WM.var(wm, n, :q_tank)

        for (i, tank) in WM.ref(wm, n, :tank)
           push!(f_2, JuMP.@NLexpression(wm.model,
               q_tank[i] * h_tank[tank["node"]]))
        end

        q_reservoir = WM.var(wm, n, :q_reservoir)

        for (i, reservoir) in WM.ref(wm, n, :reservoir)
            push!(f_2, JuMP.@NLexpression(wm.model,
                q_reservoir[i] * h_reservoir[reservoir["node"]]))
        end
    end

    JuMP.@NLobjective(wm.model, WM.JuMP.MOI.MIN_SENSE,
        sum(f_1[k] for k in 1:length(f_1)) -
        sum(f_2[k] for k in 1:length(f_2)))
end

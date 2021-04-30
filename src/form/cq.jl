"Create flow-related variables common to all directed flow models for node-connecting components."
function WM.variable_flow(wm::AbstractCQModel; nw::Int = WM._IM.nw_id_default, bounded::Bool = false, report::Bool = true)
    for name in ["des_pipe", "pipe", "pump", "regulator", "short_pipe", "valve"]
        # Create directed flow (`qp` and `qn`) variables for each component.
        WM._variable_component_flow(wm, name; nw = nw, bounded = bounded, report = report)
    end
end


function constraint_flow_conservation(
    wm::AbstractCQModel, n::Int, i::Int, pipe_fr::Array{Int64,1},
    pipe_to::Array{Int64,1}, des_pipe_fr::Array{Int64,1}, des_pipe_to::Array{Int64,1},
    pump_fr::Array{Int64,1}, pump_to::Array{Int64,1}, regulator_fr::Array{Int64,1},
    regulator_to::Array{Int64,1}, short_pipe_fr::Array{Int64,1},
    short_pipe_to::Array{Int64,1}, valve_fr::Array{Int64,1}, valve_to::Array{Int64,1},
    reservoirs::Array{Int64,1}, tanks::Array{Int64,1}, dispatchable_demands::Array{Int64,1},
    fixed_demand::Float64)
    # Collect flow variable references per component.
    q_pipe, q_des_pipe = WM.var(wm, n, :q_pipe), WM.var(wm, n, :q_des_pipe)
    q_pump, q_regulator = WM.var(wm, n, :q_pump), WM.var(wm, n, :q_regulator)
    q_short_pipe, q_valve = WM.var(wm, n, :q_short_pipe), WM.var(wm, n, :q_valve)
    q_reservoir, q_tank = WM.var(wm, n, :q_reservoir), WM.var(wm, n, :q_tank)
    q_demand = WM.var(wm, n, :q_demand)

    # Add the flow conservation constraint.
    con(wm, n, :flow_conservation)[i] = JuMP.@constraint(wm.model, -
        sum(q_pipe[a] for a in pipe_fr) + sum(q_pipe[a] for a in pipe_to) -
        sum(q_des_pipe[a] for a in des_pipe_fr) + sum(q_des_pipe[a] for a in des_pipe_to) -
        sum(z_pump[a] * q_pump[a] for a in pump_fr) +
        sum(z_pump[a] * q_pump[a] for a in pump_to) -
        sum(z_regulator[a] * q_regulator[a] for a in regulator_fr) +
        sum(z_regulator[a] * q_regulator[a] for a in regulator_to) -
        sum(q_short_pipe[a] for a in short_pipe_fr) +
        sum(q_short_pipe[a] for a in short_pipe_to) -
        sum(z_valve[a] * q_valve[a] for a in valve_fr) +
        sum(z_valve[a] * q_valve[a] for a in valve_to) == -
        sum(q_reservoir[k] for k in reservoirs) - sum(q_tank[k] for k in tanks) +
        sum(q_demand[k] for k in dispatchable_demands) + fixed_demand)
end


""
function parameter_pump_indicator(wm::AbstractCQModel; nw::Int = WM.nw_id_default)
    JuMP.@NLparameter(wm.model, z_pump[a in WM.ids(wm, nw, :pump)] == 0.0)
    wm.model[:z_pump] = z_pump
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
    wm.model[:z_regulator] = z_regulator
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
    wm.model[:z_valve] = z_valve
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
    wm.model[:h_tank] = h_tank
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
    wm.model[:h_reservoir] = h_reservoir
end


""
function update_reservoir_head(wm::AbstractCQModel, data::Dict{Int, Float64}; nw::Int = WM.nw_id_default)
    for (i, reservoir) in WM.ref(wm, nw, :reservoir)
        node = WM.ref(wm, nw, :node, reservoir["node"])
        JuMP.set_value(wm.model[:h_reservoir][node["index"]], node["head_nominal"])
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
    base_length = get(wm.data, "base_length", 1.0)
    base_time = get(wm.data, "base_time", 1.0)
    alpha = WM._get_alpha_min_1(wm) + 1.0
    
    pipe_type = wm.ref[:it][WM.wm_it_sym][:head_loss]
    viscosity = wm.ref[:it][WM.wm_it_sym][:viscosity]

    # Get NLparameter references.
    h_tank = wm.model[:h_tank]
    h_reservoir = wm.model[:h_reservoir]
    z_pump = wm.model[:z_pump]
    z_regulator = wm.model[:z_regulator]
    z_valve = wm.model[:z_valve]

    f_1 = Array{Any, 1}([0.0])
    f_2 = Array{Any, 1}([0.0])
    f_3 = Array{Any, 1}([0.0])

    for (n, network) in WM.nws(wm)
        # Get pipe flow variables.
        qp_pipe = WM.var(wm, n, :qp_pipe)
        qn_pipe = WM.var(wm, n, :qn_pipe)

        # Get design pipe flow variables.
        qp_des_pipe = WM.var(wm, n, :qp_des_pipe)
        qn_des_pipe = WM.var(wm, n, :qn_des_pipe)

        # Get valve flow variables.
        qp_valve = WM.var(wm, n, :qp_valve)
        qn_valve = WM.var(wm, n, :qn_valve)

        # Get pump flow and head difference variables.
        q_tank = WM.var(wm, n, :q_tank)
        qp_pump = WM.var(wm, n, :qp_pump)
        
        for (a, pipe) in WM.ref(wm, n, :pipe)
            L_x_r = pipe["length"] * WM._calc_pipe_resistance(
                pipe, pipe_type, viscosity, base_length, base_time)
            
            push!(f_1, JuMP.@NLexpression(wm.model, L_x_r * head_loss(qp_pipe[a])))
            push!(f_1, JuMP.@NLexpression(wm.model, L_x_r * head_loss(qn_pipe[a])))
        end

        for (a, des_pipe) in WM.ref(wm, n, :des_pipe)
            L_x_r = des_pipe["length"] * WM._calc_pipe_resistance(
                des_pipe, pipe_type, viscosity, base_length, base_time)
            
            push!(f_1, JuMP.@NLexpression(wm.model, L_x_r * head_loss(qp_des_pipe[a])))
            push!(f_1, JuMP.@NLexpression(wm.model, L_x_r * head_loss(qn_des_pipe[a])))
        end
        
        for (a, pump) in WM.ref(wm, n, :pump)
            @assert pump["head_curve_form"] in [WM.PUMP_QUADRATIC,
                WM.PUMP_BEST_EFFICIENCY_POINT, WM.PUMP_LINEAR_POWER]

            c = WM._calc_head_curve_coefficients(pump)
            
            push!(f_2, JuMP.@NLexpression(wm.model,
                z_pump[a] * c[1] / 3.0 * qp_pump[a]^3 +
                z_pump[a] * 0.5 * c[2] * qp_pump[a]^2 +
                z_pump[a] * c[3] * qp_pump[a]))
        end

        for (i, res) in WM.ref(wm, n, :reservoir)
            head = WM.ref(wm, n, :node, res["node"])["head_nominal"]

            for comp in [:des_pipe, :pipe, :pump, :regulator, :short_pipe, :valve]
                qp = WM.var(wm, n, Symbol("qp_" * string(comp)))
                qn = WM.var(wm, n, Symbol("qn_" * string(comp)))

                for a in WM.ref(wm, n, Symbol(string(comp) * "_fr"), res["node"])
                    reservoir_head_times_flow = JuMP.@NLexpression(
                        wm.model, h_reservoir[res["node"]] * (qp[a] - qn[a]))
                    push!(f_2, JuMP.@NLexpression(wm.model, reservoir_head_times_flow))
                end

                for a in WM.ref(wm, n, Symbol(string(comp) * "_to"), res["node"])
                    reservoir_head_times_flow = JuMP.@NLexpression(
                        wm.model, h_reservoir[res["node"]] * (qn[a] - qp[a]))
                    push!(f_2, JuMP.@NLexpression(wm.model, reservoir_head_times_flow))
                end
            end
        end

        for (i, tank) in WM.ref(wm, n, :tank)
            push!(f_3, JuMP.@NLexpression(wm.model,
                -q_tank[i] * h_tank[tank["node"]]))
        end
    end

    JuMP.@NLobjective(wm.model, WM._MOI.MIN_SENSE,
        sum(f_1[k] for k in 1:length(f_1)) -
        sum(f_2[k] for k in 1:length(f_2)) +
        sum(f_3[k] for k in 1:length(f_3)))
end

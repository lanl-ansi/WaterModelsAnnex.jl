function add_price_cuts!(wm::WM.AbstractWaterModel)
    nw_ids = sort(collect(WM.nw_ids(wm)))[1:end-1]

    for partition in collect(Iterators.partition(nw_ids, 2))
        for pipe_id in WM.ids(wm, 1, :pipe)
            q_pipe = WM.var.(Ref(wm), partition, :q_pipe, pipe_id)
            JuMP.@constraint(wm.model, q_pipe[1] == q_pipe[2])

            y_pipe = WM.var.(Ref(wm), partition, :y_pipe, pipe_id)
            JuMP.@constraint(wm.model, y_pipe[1] == y_pipe[2])
        end

        for short_pipe_id in WM.ids(wm, 1, :short_pipe)
            q_short_pipe = WM.var.(Ref(wm), partition, :q_short_pipe, short_pipe_id)
            JuMP.@constraint(wm.model, q_short_pipe[1] == q_short_pipe[2])

            y_short_pipe = WM.var.(Ref(wm), partition, :y_short_pipe, short_pipe_id)
            JuMP.@constraint(wm.model, y_short_pipe[1] == y_short_pipe[2])
        end

        for pump_id in WM.ids(wm, 1, :pump)
            z_pump = WM.var.(Ref(wm), partition, :z_pump, pump_id)
            JuMP.@constraint(wm.model, z_pump[1] == z_pump[2])

            y_pump = WM.var.(Ref(wm), partition, :y_pump, pump_id)
            JuMP.@constraint(wm.model, y_pump[1] == y_pump[2])

            q_pump = WM.var.(Ref(wm), partition, :q_pump, pump_id)
            JuMP.@constraint(wm.model, q_pump[1] == q_pump[2])
        end

        for valve_id in WM.ids(wm, 1, :valve)
            y_valve = WM.var.(Ref(wm), partition, :y_valve, valve_id)
            JuMP.@constraint(wm.model, y_valve[1] == y_valve[2])

            z_valve = WM.var.(Ref(wm), partition, :z_valve, valve_id)
            JuMP.@constraint(wm.model, z_valve[1] == z_valve[2])

            q_valve = WM.var.(Ref(wm), partition, :q_valve, valve_id)
            JuMP.@constraint(wm.model, q_valve[1] == q_valve[2])
        end
    end

    for nw in nw_ids
        for (a, pump) in WM.ref(wm, nw, :pump)
            z_pump_current = WM.var(wm, nw, :z_pump)
            energy_price = WM.ref(wm, nw, :pump, a, "energy_price")
            energy_prices = WM.ref.(Ref(wm), nw_ids, :pump, a, "energy_price")
            nw_ids_price = nw_ids[findall(x -> x < energy_price, energy_prices)]
            
            # if length(nw_ids_price) > 0
            #     z_pump_sum = sum(WM.var(wm, n, :z_pump, a) for n in nw_ids_price)
            #     JuMP.@constraint(wm.model, z_pump_sum >= z_pump_current)
            # end

            for nw_id_price in setdiff(Set(nw_ids_price), [nw])
                # println(energy_prices[nw_id_price], " ", energy_price)
                z_pump_nw_price = WM.var(wm, nw_id_price, :z_pump)
                JuMP.@constraint(wm.model, sum(z_pump_nw_price) >= sum(z_pump_current))
            end
        end
    end
end
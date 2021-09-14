function solve_obbt_owf_volume!(
    data::Dict{String,<:Any}, optimizer; use_relaxed_network::Bool = true,
    model_type::Type = WM.PWLRDWaterModel, time_limit::Float64 = 3600.0,
    upper_bound::Float64 = Inf, upper_bound_constraint::Bool = false,
    max_iter::Int = 100, solve_relaxed::Bool = true, limit_problems::Bool = false,
    error_tolerance::Float64 = 1.0, flow_tolerance::Float64 = 1.0e-4, kwargs...)
    # Print a message with relevant algorithm limit information.
    message = "[OBBT] Maximum time limit set to default value of $(time_limit) seconds."
    WM.Memento.info(WM._LOGGER, message)

    # Relax the network (e.g., make nodal components dispatchable) if requested.
    use_relaxed_network && WM.relax_network!(data)

    # Set the problem specification that will be used for bound tightening.
    build_type = WM._IM.ismultinetwork(WM.get_wm_data(data)) ? build_mn_owf_volume : build_wf

    # Check for keyword argument inconsistencies.
    WM._check_obbt_options(upper_bound, upper_bound_constraint)

    # Instantiate the bound tightening model and relax integrality, if requested.
    wms = [WM.instantiate_model(data, model_type, build_type) for i in 1:Threads.nthreads()]

    if upper_bound_constraint
        map(x -> WM._constraint_obj_bound(x, upper_bound), wms)
    end

    if solve_relaxed
        # Relax the binary variables if requested.
        map(x -> WM.relax_all_binary_variables!(x), wms)
    end

    # Set the optimizer for the bound tightening model.
    map(x -> JuMP.set_optimizer(x.model, optimizer), wms)

    # Collect all problems.
    bound_problems = _get_bound_problems(wms[1]; limit = limit_problems)
    WM._update_data_bounds!(data, bound_problems) # Populate data with bounds.

    # Log widths.
    bound_width_msg = WM._log_bound_widths(data)
    WM.Memento.info(WM._LOGGER, "[OBBT] Initial bound widths: $(bound_width_msg).")
    terminate, time_elapsed = false, 0.0

    # Set up algorithm metadata.
    current_iteration = 1
    terminate = current_iteration >= max_iter

    while any([x.changed for x in bound_problems]) && !terminate
        # Obtain new candidate bounds, update bounds, and update the data.
        vals = zeros(length(bound_problems))

        time_elapsed += @elapsed Threads.@threads for i in 1:length(bound_problems)
            vals[i] = WM._solve_bound_problem!(wms[Threads.threadid()], bound_problems[i])
            WM._set_new_bound!(bound_problems[i], vals[i])
        end

        time_elapsed > time_limit && ((terminate = true) && break)
        WM._update_data_bounds!(data, bound_problems)
        WM.set_flow_partitions_si!(data, error_tolerance, flow_tolerance)
        !terminate && WM._clean_bound_problems!(bound_problems, vals)

        # Log widths.
        bound_width_msg = WM._log_bound_widths(data)
        message = "[OBBT] Iteration $(current_iteration) bound widths: $(bound_width_msg)."
        WM.Memento.info(WM._LOGGER, message)

        # Update algorithm metadata.
        current_iteration += 1

        # Set up the next optimization problem using the new bounds.
        wms = [WM.instantiate_model(data, model_type, build_type) for i in 1:Threads.nthreads()]

        if upper_bound_constraint
            map(x -> WM._constraint_obj_bound(x, upper_bound), wms)
        end

        if solve_relaxed
            # Relax the binary variables if requested.
            map(x -> WM.relax_all_binary_variables!(x), wms)
        end

        # Set the optimizer for the bound tightening model.
        map(x -> JuMP.set_optimizer(x.model, optimizer), wms)

        # Set the termination variable if max iterations is exceeded.
        current_iteration >= max_iter && (terminate = true)
    end

    if use_relaxed_network
        WM._fix_demands!(data)
        WM._fix_tanks!(data)
        WM._fix_reservoirs!(data)
    end

    time_elapsed_rounded = round(time_elapsed; digits = 2)
    WM.Memento.info(WM._LOGGER, "[OBBT] Completed in $(time_elapsed_rounded) seconds.")
end


function solve_obbt_owf_switching!(
    data::Dict{String,<:Any}, optimizer; use_relaxed_network::Bool = true,
    model_type::Type = WM.PWLRDWaterModel, time_limit::Float64 = 3600.0,
    upper_bound::Float64 = Inf, upper_bound_constraint::Bool = false,
    max_iter::Int = 100, solve_relaxed::Bool = true, cuts::Vector{WM._PairwiseCut} = [],
    flow_partition_func::Function = x -> WM.set_flow_partitions_si!(x, 1.0, 1.0e-4),
    limit_problems::Bool = false, kwargs...)
    # Print a message with relevant algorithm limit information.
    message = "[OBBT] Maximum time limit set to default value of $(time_limit) seconds."
    WM.Memento.info(WM._LOGGER, message)

    # Relax the network (e.g., make nodal components dispatchable) if requested.
    use_relaxed_network && WM.relax_network!(data)

    # Set the problem specification that will be used for bound tightening.
    build_type = WM._IM.ismultinetwork(WM.get_wm_data(data)) ? WM.build_mn_owf_switching : WM.build_wf

    # Check for keyword argument inconsistencies.
    WM._check_obbt_options(upper_bound, upper_bound_constraint)

    # Instantiate the bound tightening model and relax integrality, if requested.
    wms = Vector{WM.AbstractWaterModel}(undef, Threads.nthreads())

    # Update WaterModels objects in parallel.
    Threads.@threads for i in 1:Threads.nthreads()
        wms[i] = WM.instantiate_model(data, model_type, build_type)

        # Add pairwise cuts to the models.
        add_pairwise_cuts(wms[i], cuts)

        if upper_bound_constraint
            WM._constraint_obj_bound(wms[i], upper_bound)
        end

        if solve_relaxed
            # Relax the binary variables if requested.
            WM.relax_all_binary_variables!(wms[i])
        end

        # Set the optimizer for the bound tightening model.
        JuMP.set_optimizer(wms[i].model, optimizer)
    end

    # Collect all problems.
    bound_problems = WM._get_bound_problems(wms[1]; limit = limit_problems)
    WM._update_data_bounds!(data, bound_problems) # Populate data with bounds.

    # Log widths.
    bound_width_msg = WM._log_bound_widths(data)
    WM.Memento.info(WM._LOGGER, "[OBBT] Initial bound widths: $(bound_width_msg).")
    terminate, time_elapsed = false, 0.0

    # Set up algorithm metadata.
    current_iteration = 1
    terminate = current_iteration >= max_iter

    while any([x.changed for x in bound_problems]) && !terminate
        # Obtain new candidate bounds, update bounds, and update the data.
        vals = zeros(length(bound_problems))

        time_elapsed += @elapsed Threads.@threads for i in 1:length(bound_problems)
            vals[i] = WM._solve_bound_problem!(wms[Threads.threadid()], bound_problems[i])
            WM._set_new_bound!(bound_problems[i], vals[i])
        end

        time_elapsed > time_limit && ((terminate = true) && break)
        WM._update_data_bounds!(data, bound_problems)
        flow_partition_func(data)
        !terminate && WM._clean_bound_problems!(bound_problems, vals)

        # Log widths.
        bound_width_msg = WM._log_bound_widths(data)
        message = "[OBBT] Iteration $(current_iteration) bound widths: $(bound_width_msg)."
        WM.Memento.info(WM._LOGGER, message)

        # Update algorithm metadata.
        current_iteration += 1

        # Update WaterModels objects in parallel.
        Threads.@threads for i in 1:Threads.nthreads()
            # Set up the next optimization problem using the new bounds.
            wms[i] = WM.instantiate_model(data, model_type, build_type)

            # Add pairwise cuts to the models.
            add_pairwise_cuts(wms[i], cuts)

            if upper_bound_constraint
                WM._constraint_obj_bound(wms[i], upper_bound)
            end

            if solve_relaxed
                # Relax the binary variables if requested.
                WM.relax_all_binary_variables!(wms[i])
            end

            # Set the optimizer for the bound tightening model.
            JuMP.set_optimizer(wms[i].model, optimizer)
        end

        # Set the termination variable if max iterations is exceeded.
        current_iteration >= max_iter && (terminate = true)
    end

    if use_relaxed_network
        WM._fix_demands!(data)
        WM._fix_tanks!(data)
        WM._fix_reservoirs!(data)
    end

    time_elapsed_rounded = round(time_elapsed; digits = 2)
    WM.Memento.info(WM._LOGGER, "[OBBT] Completed in $(time_elapsed_rounded) seconds.")
end

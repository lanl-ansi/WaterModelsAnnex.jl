function get_owf_heuristic_callback(wm::WM.AbstractWaterModel)
    return function callback_function(cb_data)
    end
end


function add_owf_heuristic_callback!(wm::WM.AbstractWaterModel)
    callback_function = get_owf_heuristic_callback(wm)
    WM._MOI.set(wm.model, WM._MOI.HeuristicCallback(), callback_function)
end
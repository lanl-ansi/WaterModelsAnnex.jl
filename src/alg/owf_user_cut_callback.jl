function get_owf_user_cut_callback(wm::WM.AbstractWaterModel)
    return function callback_function(cb_data)
    end
end


function add_owf_user_cut_callback!(wm::WM.AbstractWaterModel)
    callback_function = get_owf_user_cut_callback(wm)
    WM._MOI.set(wm.model, WM._MOI.UserCutCallback(), callback_function)
end
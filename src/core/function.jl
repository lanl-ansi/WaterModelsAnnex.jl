#################################################################################
# This file defines the nonlinear head loss functions for water systems models. #
#################################################################################


function _fq_cd(alpha::Float64)
    return x -> 1.0 / (2.0 + alpha) * (x*x)^(1.0 + 0.5*alpha)
end


function _dfq_cd(alpha::Float64)
    return x -> (x*x)^(0.5 + 0.5*alpha)
end


function _d2fq_cd(alpha::Float64)
    return x -> (1.0 + alpha) * sign(x) * (x*x)^(0.5*alpha)
end


function _fdh_cd(alpha::Float64)
    return x -> (alpha + 1.0) / (alpha + 2.0) * (x * x)^(0.5 + 0.5 / (1.0 + alpha))
end


function _dfdh_cd(alpha::Float64)
    return x -> (x*x)^(0.5 / (alpha + 1.0))
end


function _d2fdh_cd(alpha::Float64)
    return x -> sign(x) * ((x*x)^(-(0.5 * alpha) / (1.0 + alpha))) / (1.0 + alpha)
end


function head_loss_args(wm::AbstractCDModel)
    alpha_m1 = WM._get_alpha_min_1(wm)
    return (:head_loss, 1, _fq_cd(alpha_m1), _dfq_cd(alpha_m1), _d2fq_cd(alpha_m1))
end


function head_loss_dh_args(wm::AbstractCDModel)
    alpha_m1 = WM._get_alpha_min_1(wm)
    return (:head_loss_dh, 1, _fdh_cd(alpha_m1), _dfdh_cd(alpha_m1), _d2fdh_cd(alpha_m1))
end


function WM._function_head_loss(wm::AbstractCDModel)
    JuMP.register(wm.model, head_loss_args(wm)...)
    JuMP.register(wm.model, head_loss_dh_args(wm)...)
end
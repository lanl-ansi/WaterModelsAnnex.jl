"Directed nonlinear, convex models."
abstract type AbstractCDModel <: WM.AbstractNCDModel end
mutable struct CDWaterModel <: AbstractCDModel WM.@wm_fields end
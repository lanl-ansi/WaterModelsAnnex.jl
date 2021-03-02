"Directed nonlinear, convex models."
abstract type AbstractCDModel <: WM.AbstractNCDModel end
mutable struct CDWaterModel <: AbstractCDModel WM.@wm_fields end

"Extended directed nonlinear, convex models."
abstract type AbstractCDXModel <: WM.AbstractCRDModel end
mutable struct CDXWaterModel <: AbstractCDXModel WM.@wm_fields end

"Extended directed linear, relaxation-based models."
abstract type AbstractLRDXModel <: WM.AbstractLRDModel end
mutable struct LRDXWaterModel <: AbstractLRDXModel WM.@wm_fields end
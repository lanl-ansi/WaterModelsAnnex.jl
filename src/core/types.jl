"Directed nonlinear, convex models."
abstract type AbstractCDModel <: WM.AbstractNCDModel end
mutable struct CDWaterModel <: AbstractCDModel WM.@wm_fields end

"Nonlinear, nonconvex models."
abstract type AbstractNCZModel <: WM.AbstractNCModel end
mutable struct NCZWaterModel <: AbstractNCZModel WM.@wm_fields end

"Directed nonlinear, convex models."
abstract type AbstractCQModel <: WM.AbstractNCDModel end
mutable struct CQWaterModel <: AbstractCQModel WM.@wm_fields end

"Extended directed nonlinear, convex models."
abstract type AbstractCDXModel <: WM.AbstractCRDModel end
mutable struct CDXWaterModel <: AbstractCDXModel WM.@wm_fields end

"Extended directed linear, relaxation-based models."
abstract type AbstractPWLRDXModel <: WM.AbstractPWLRDModel end
mutable struct PWLRDXWaterModel <: AbstractPWLRDXModel WM.@wm_fields end
using WaterModelsAnnex
import InfrastructureModels
import Memento
import WaterModels

import Cbc
import Ipopt
import JuMP
import JSON

import MathOptInterface
const MOI = MathOptInterface
const MOIU = MathOptInterface.Utilities
const WMs = WaterModels

using Test

# Suppress warnings during testing.
Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
Memento.setlevel!(Memento.getlogger(WaterModels), "error")
Memento.setlevel!(Memento.getlogger(WaterModelsAnnex), "error")

# Default MIP and NLP optimizers.
const cbc = JuMP.with_optimizer(Cbc.Optimizer, logLevel=0)
const ipopt = JuMP.with_optimizer(Ipopt.Optimizer, tol=1.0e-9, acceptable_tol=1.0e-9, max_iter=9999, print_level=0)

@testset "WaterModelsAnnex" begin

    include("alg.jl")

end

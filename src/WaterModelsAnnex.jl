module WaterModelsAnnex

import BlockDecomposition
import Coluna
import Gurobi
import JuMP
import LinearAlgebra
import PolyhedralRelaxations
import WaterModels
import LightGraphs
import MetaGraphs

const BD = BlockDecomposition
const WM = WaterModels
const JuMP = WM.JuMP
const MOI = WM.MathOptInterface
const MOIU = WM.MathOptInterface.Utilities
const LOGGER = WM.Memento.getlogger(WM)

# Register the module level logger at runtime so that it can be accessed via
# `getlogger(WaterModelsAnnex)` NOTE: If this line is not included then the
# precompiled `WaterModelsAnnex.LOGGER` won't be registered at runtime.
__init__() = WM.Memento.register(LOGGER)

"Suppresses information and warning messages output by WaterModels. For fine-grained control use the Memento package."
function silence()
    WM.Memento.info(LOGGER, "Suppressing information and warning messages for the rest of this session  Use the Memento package for fine-grained control of logging.")
    WM.Memento.setlevel!(Memento.getlogger(WM._IM), "error")
    WM.Memento.setlevel!(Memento.getlogger(WM), "error")
end

include("core/types.jl")
include("core/function.jl")
include("core/constraint.jl")
include("core/constraint_template.jl")
include("core/control_setting.jl")
include("core/simulation_result.jl")
include("core/variable.jl")

include("form/cd.jl")
include("form/cq.jl")
include("form/cdx.jl")
include("form/lrdx.jl")

include("prob/wf.jl")
include("prob/owf.jl")
include("prob/owfh.jl")

include("alg/set_breakpoints.jl")
include("alg/dantzig_wolfe.jl")
include("alg/simulate.jl")
include("alg/heuristic_relaxed_milp.jl")
include("alg/owf_lazy_cut_callback.jl")
include("alg/owf_user_cut_callback.jl")
include("alg/owf_heuristic_callback.jl")

include("alg/heuristic_cuts.jl")
include("alg/compute_source_pumps.jl")
include("alg/pump_volume_cuts.jl")
include("alg/solve_owf.jl")
include("alg/heuristic_linear_program.jl")
include("alg/heuristic_master_program.jl")
include("alg/repair_schedule.jl")

include("util/bound_problem.jl")
include("util/graph.jl")
include("util/create_schedules.jl")
include("util/warm_start.jl")
include("util/obbt.jl")
include("core/export.jl")

end

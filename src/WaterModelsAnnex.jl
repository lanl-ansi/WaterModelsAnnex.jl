module WaterModelsAnnex

import BlockDecomposition
import Coluna
import JuMP
import WaterModels
import LightGraphs

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
include("core/constraint_template.jl")

include("form/cd.jl")
include("form/cdx.jl")
include("form/lrdx.jl")

include("prob/wf.jl")
include("prob/owfh.jl")

include("core/export.jl")

include("alg/dantzig_wolfe.jl")
include("alg/simulate.jl")
include("alg/owf_lazy_cut_callback.jl")
include("alg/owf_user_cut_callback.jl")
include("alg/owf_heuristic_callback.jl")

include("alg/solve_owf.jl")

include("util/create_schedules.jl")
include("core/export.jl")

end

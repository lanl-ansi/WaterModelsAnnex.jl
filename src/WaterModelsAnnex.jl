module WaterModelsAnnex

using WaterModels

const WM = WaterModels
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

include("alg/solve-owf.jl")

#include("alg/ne-add_tasseff.jl")
#include("alg/ne-global.jl")
#include("alg/ne-decomp_cut.jl")
#include("alg/ne-humpola_cut.jl")
#include("alg/ne-no_good_cut.jl")
#include("alg/ne-tasseff_cut.jl")
#include("util/check_physical_solution.jl")
#include("util/get_physical_solution.jl")

end

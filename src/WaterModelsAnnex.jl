module WaterModelsAnnex

import InfrastructureModels
import JSON
import JuMP
import Memento
using WaterModels

const LOGGER = Memento.getlogger(WaterModels)

const WMs = WaterModels

# Register the module level logger at runtime so that it can be accessed via
# `getlogger(WaterModelsAnnex)` NOTE: If this line is not included then the
# precompiled `WaterModelsAnnex.LOGGER` won't be registered at runtime.
__init__() = Memento.register(LOGGER)

"Suppresses information and warning messages output by WaterModels. For fine-grained control use the Memento package."
function silence()
    Memento.info(LOGGER, "Suppressing information and warning messages for the rest of this session  Use the Memento package for fine-grained control of logging.")
    Memento.setlevel!(Memento.getlogger(InfrastructureModels), "error")
    Memento.setlevel!(Memento.getlogger(WaterModels), "error")
end

include("alg/ne-no_good.jl")
#include("alg/ne-repair_solution.jl")
#include("alg/ne-find_initial_solution.jl")
#include("alg/ne-benders.jl")
#include("alg/ne-raghunathan.jl")
#include("alg/ne-tasseff.jl")

end

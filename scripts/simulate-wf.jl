import Gurobi
import Ipopt
import JuMP

using Revise
using WaterModelsAnnex

const WM = WaterModelsAnnex.WM

# Specify the network and multinetwork datasets.
network_path = "WaterModels.jl/examples/data/epanet/van_zyl.inp";
network = WM.parse_file(network_path);

# Set up Gurobi solvers.
env = Gurobi.Env();
gurobi_1 = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "MIPGap" => 0.0);
gurobi_2 = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "TimeLimit" => 20.0);

# Tighten bounds of variables in the network.
WM.solve_obbt_owf!(network, gurobi_1; model_type = WM.PWLRDWaterModel,
    solve_relaxed = false, time_limit = 3600.0, max_iter = 100,
    ext = Dict(:pipe_breakpoints => 10, :pump_breakpoints => 10));

# Create the multinetwork version of the network.
network_mn = WM.make_multinetwork(network);

# Solve a relaxation of the OWF problem.
ext = Dict{Symbol, Any}(:pipe_breakpoints => 10, :pump_breakpoints => 10);
#result = WM.solve_mn_owf(network_mn, WM.PWLRDWaterModel, gurobi_2; ext = ext);
result = solve_mn_owfh(network_mn, WM.PWLRDWaterModel, gurobi_2; ext = ext);

ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "acceptable_tol" => 1.0e-8, "sb" => "yes");
result = WaterModelsAnnex.simulate!(network, result, ipopt);
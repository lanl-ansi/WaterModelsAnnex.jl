import Gurobi
import Ipopt
import JuMP

using Revise
using WaterModelsAnnex

const WM = WaterModelsAnnex.WM

# Specify the network and multinetwork datasets.
network_path = "examples/data/epanet/van_zyl.inp";
network = WM.parse_file(network_path);

# Set up Gurobi solvers.
env = Gurobi.Env();
gurobi = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "MIPGap" => 0.0);
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes");

schedules = calc_possible_schedules(network, ipopt);
weights = solve_heuristic_problem(network, schedules, gurobi);
sched_sol = solve_heuristic_master(network, schedules, weights, ipopt, gurobi);

# network = WM.parse_file(network_path);
# network_mn = WM.make_multinetwork(network);
# ext = Dict{Symbol, Any}(:pipe_breakpoints => 10, :pump_breakpoints => 10);
# result_rel = WM.solve_mn_owf(network_mn, WM.PWLRDWaterModel, gurobi_2; relax_integrality = true);

# # Reload the network data.
# network = WM.parse_file(network_path);
# result = WaterModelsAnnex.simulate!(network, schedules, sched_sol, result_rel, ipopt);

# # Tighten bounds of variables in the network.
# WM.solve_obbt_owf!(network, gurobi_1; model_type = WM.PWLRDWaterModel,
#     solve_relaxed = false, time_limit = 3600.0, max_iter = 100,
#     ext = Dict(:pipe_breakpoints => 10, :pump_breakpoints => 10));

# # Solve a relaxation of the OWF problem.
# ext = Dict{Symbol, Any}(:pipe_breakpoints => 10, :pump_breakpoints => 10);
# #result = WM.solve_mn_owf(network_mn, WM.PWLRDWaterModel, gurobi_2; ext = ext);
# result = solve_mn_owfh(network_mn, WM.PWLRDWaterModel, gurobi_2; ext = ext);

# ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "acceptable_tol" => 1.0e-8, "sb" => "yes");
# result = WaterModelsAnnex.simulate!(network, result, ipopt);
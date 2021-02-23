import Gurobi
import Ipopt
import JuMP
import Memento
import Logging

using Revise
using WaterModelsAnnex

const WM = WaterModelsAnnex.WM
Memento.setlevel!(Memento.getlogger(WM._IM), "error")

# Set up Gurobi solvers.
env = Gurobi.Env();
gurobi = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "MIPGap" => 0.0);
gurobi_owf = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "MIPGap" => 0.0, "TimeLimit" => 120.0);
gurobi_short = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "MIPGap" => 0.0, "TimeLimit" => 5.0);
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes", "linear_solver" => "ma57");

# Suppress warnings during testing.
Memento.setlevel!(Memento.getlogger(WM), "info");
Memento.setlevel!(Memento.getlogger(WM._IM), "error");
Logging.disable_logging(Logging.Info);

# Specify the network and multinetwork datasets.
network_path = "examples/data/epanet/van_zyl.inp"
network = WM.parse_file(network_path; skip_correct = true);
modifications = WM.parse_file("modifications.json"; skip_correct = true);
WM._IM.update_data!(network, modifications);
WM.correct_network_data!(network);

schedules = calc_possible_schedules(network, ipopt);
weights = solve_heuristic_problem(network, schedules, gurobi);
result_heur = solve_heuristic_master(network, schedules, weights, ipopt, gurobi);

network = WM.parse_file(network_path; skip_correct = true);
modifications = WM.parse_file("modifications.json"; skip_correct = true);
WM._IM.update_data!(network, modifications);
WM.correct_network_data!(network);

# Specify the network and multinetwork datasets.
WM._IM.update_data!(network_mn, result_heur["solution"]);


network_path = "examples/data/epanet/van_zyl.inp"
network = WM.parse_file(network_path);

# Tighten bounds of variables in the network.
WM.solve_obbt_owf!(network, gurobi; model_type = WM.PWLRDWaterModel,
    solve_relaxed = false, time_limit = 3600.0, max_iter = 100,
    ext = Dict(:pipe_breakpoints => 10, :pump_breakpoints => 10));

network_mn = WM.make_multinetwork(network);
ext = Dict{Symbol, Any}(:pipe_breakpoints => 6, :pump_breakpoints => 6);
result = WM.solve_mn_owf(network_mn, WM.PWLRDWaterModel, gurobi_owf; ext = ext);

WaterModelsAnnex._update_tank_time_series!(network, result)
schedules = calc_possible_schedules(network, ipopt);
weights = solve_heuristic_problem(network, schedules, gurobi);
result_heur = solve_heuristic_master(network, schedules, weights, ipopt, gurobi);

WM._IM.update_data!(network_mn, result_heur["solution"]);
WM.set_start_all!(network_mn);
# WM.turn_on_all_components!(network_mn);
result = WM.solve_mn_owf(network_mn, WM.LRDWaterModel, gurobi_short; ext = ext);

# # Reload the network data.
# network = WM.parse_file(network_path);
# result = WaterModelsAnnex.simulate!(network, schedules, sched_sol, result_rel, ipopt);

# # Solve a relaxation of the OWF problem.
# ext = Dict{Symbol, Any}(:pipe_breakpoints => 10, :pump_breakpoints => 10);
# #result = WM.solve_mn_owf(network_mn, WM.PWLRDWaterModel, gurobi_2; ext = ext);
# result = solve_mn_owfh(network_mn, WM.PWLRDWaterModel, gurobi_2; ext = ext);

# ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "acceptable_tol" => 1.0e-8, "sb" => "yes");
# result = WaterModelsAnnex.simulate!(network, result, ipopt);
import Gurobi
import Ipopt
import JuMP
import Memento
import Logging
import Plots

using Revise
using WaterModelsAnnex

Memento.setlevel!(Memento.getlogger(WaterModelsAnnex.WM._IM), "error")

# Set up optimization solvers.
env = Gurobi.Env();
gurobi = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "OutputFlag" => 0, "MIPGap" => 0.0);
gurobi_owf = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "MIPGap" => 0.0, "TimeLimit" => 60.0);
gurobi_short = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "MIPGap" => 0.0, "TimeLimit" => 5.0);
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes");

# Parse example network data.
network_path = "examples/data/epanet/van_zyl.inp"
network = WM.parse_file(network_path);

# Tighten bounds of variables in the network.
WM.solve_obbt_owf!(network, gurobi; model_type = WM.PWLRDWaterModel,
    solve_relaxed = false, time_limit = 3600.0, max_iter = 100,
    ext = Dict(:pipe_breakpoints => 10, :pump_breakpoints => 10));

# Solve a relaxed version of the problem.
network_milp = deepcopy(network)
network_mn = WM.make_multinetwork(network_milp);
ext = Dict{Symbol, Any}(:pipe_breakpoints => 10, :pump_breakpoints => 10);
result = WM.solve_mn_owf(network_mn, WM.PWLRDWaterModel, gurobi_owf; ext = ext, relax_integrality = true);
result_milp = deepcopy(result)
result_milp_sim = WaterModelsAnnex.simulate!(network_milp, result_milp, ipopt)

# Update the time series for tank levels.
network_heur = deepcopy(network)
WaterModelsAnnex._update_tank_time_series!(network_heur, result);
schedules = calc_possible_schedules(network_heur, ipopt);
weights = solve_heuristic_problem(network_heur, schedules, ipopt);
result_heur = solve_heuristic_master(network_heur, schedules, weights, ipopt, gurobi);
result_heur_sim = WaterModelsAnnex.simulate!(network_heur, result_heur, ipopt)

wm = WM.instantiate_model(network_mn, CDWaterModel, WM.build_mn_wf);
schedules = create_all_schedules(wm);
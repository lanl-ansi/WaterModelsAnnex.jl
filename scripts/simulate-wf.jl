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
gurobi_owf = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "MIPGap" => 0.0, "TimeLimit" => 60.0);
gurobi_short = JuMP.optimizer_with_attributes(() -> Gurobi.Optimizer(env), "MIPGap" => 0.0, "TimeLimit" => 5.0);
ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "sb" => "yes", "linear_solver" => "mumps");

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
ext = Dict{Symbol, Any}(:pipe_breakpoints => 3, :pump_breakpoints => 3);
result = WM.solve_mn_owf(network_mn, WM.PWLRDWaterModel, gurobi_owf; ext = ext);

result_copy = deepcopy(result);
result_copy_milp = deepcopy(result);

network_copy = deepcopy(network);
result_milp_sim = simulate!(network_copy, result_copy_milp, ipopt);

WaterModelsAnnex._update_tank_time_series!(network_copy, result_copy)
schedules = calc_possible_schedules(network_copy, ipopt);
weights = solve_heuristic_problem(network_copy, schedules, gurobi);
result_heur = solve_heuristic_master(network_copy, schedules, weights, ipopt, gurobi);

result_fixed_sim = simulate!(network_copy, result_heur, ipopt);

tank_1_max_vol = ones(24) * 0.25 * pi * network_copy["tank"]["1"]["diameter"]^2 * network_copy["tank"]["1"]["max_level"]
tank_1_milp = [result["solution"]["nw"][string(n)]["tank"]["1"]["V"] for n in 1:24];
tank_1_milp_sim = [result_milp_sim["solution"]["nw"][string(n)]["tank"]["1"]["V"] for n in 1:24];
tank_1_fixed_sim = [result_fixed_sim["solution"]["nw"][string(n)]["tank"]["1"]["V"] for n in 1:24];
tank_1_plot_milp = plot(1:4, [tank_1_milp[1:4], tank_1_milp_sim[1:4], tank_1_max_vol[1:4]], xlabel = "Time (h)", ylabel = "Tank 1 Volume", markershape = :circle, markersize = 3, label = ["MILP" "MILP Sim" "Max Volume"], legend=(0.2, 0.4));
tank_1_plot = plot(1:24, [tank_1_milp, tank_1_fixed_sim], xlabel = "Time (h)", ylabel = "Tank 1 Volume", markershape = :circle, markersize = 3, label = ["MILP" "Repaired"], legendfontsize = 8, legend=(0.15, 0.35));

tank_2_max_vol = ones(24) * 0.25 * pi * network_copy["tank"]["2"]["diameter"]^2 * network_copy["tank"]["2"]["max_level"]
tank_2_milp = [result["solution"]["nw"][string(n)]["tank"]["2"]["V"] for n in 1:24];
tank_2_milp_sim = [result_milp_sim["solution"]["nw"][string(n)]["tank"]["2"]["V"] for n in 1:24];
tank_2_fixed_sim = [result_fixed_sim["solution"]["nw"][string(n)]["tank"]["2"]["V"] for n in 1:24];
tank_2_plot_milp = plot(1:4, [tank_2_milp[1:4], tank_2_milp_sim[1:4], tank_2_max_vol[1:4]], xlabel = "Time (h)", ylabel = "Tank 2 Volume", markershape = :circle, markersize = 3, legend = false);
tank_2_plot = plot(1:24, [tank_2_milp, tank_2_fixed_sim], xlabel = "Time (h)", ylabel = "Tank 2 Volume", markershape = :circle, markersize = 3, legend = false);

pump_1_milp = [result["solution"]["nw"][string(n)]["pump"]["1"]["status"] for n in 1:24];
pump_1_fixed_sim = [result_fixed_sim["solution"]["nw"][string(n)]["pump"]["1"]["status"] for n in 1:24];
pump_1_plot = plot(1:24, [pump_1_milp, pump_1_fixed_sim], xlabel = "Time (h)", ylabel = "Pump 1 Status", markershape = :circle, markersize = 3, label = ["MILP" "Repaired"], legendfontsize = 8, legend = (0.83, 0.5));

pump_2_milp = [result["solution"]["nw"][string(n)]["pump"]["2"]["status"] for n in 1:24];
pump_2_fixed_sim = [result_fixed_sim["solution"]["nw"][string(n)]["pump"]["2"]["status"] for n in 1:24];
pump_2_plot = plot(1:24, [pump_2_milp, pump_2_fixed_sim], xlabel = "Time (h)", ylabel = "Pump 2 Status", markershape = :circle, markersize = 3, legend = false);

pump_3_milp = [result["solution"]["nw"][string(n)]["pump"]["3"]["status"] for n in 1:24];
pump_3_fixed_sim = [result_fixed_sim["solution"]["nw"][string(n)]["pump"]["3"]["status"] for n in 1:24];
pump_3_plot = plot(1:24, [pump_3_milp, pump_3_fixed_sim], xlabel = "Time (h)", ylabel = "Pump 3 Status", markershape = :circle, markersize = 3, legend = false);

plot(pump_1_plot, pump_2_plot, pump_3_plot, layout = (3, 1), legendfontsize = 4)

valve_15_milp = [result["solution"]["nw"][string(n)]["valve"]["15"]["status"] for n in 1:24];
valve_15_fixed_sim = [result_fixed_sim["solution"]["nw"][string(n)]["valve"]["15"]["status"] for n in 1:24];
plot(1:24, [valve_15_milp, valve_15_fixed_sim], xlabel = "Time (h)", ylabel = "Valve 15 Status", markershape = :circle, markersize = 3, label = ["MILP" "Repaired"])

plot(tank_1_plot_milp, tank_2_plot_milp, layout = (2, 1))
plot(pump_1_plot, pump_2_plot, pump_3_plot, layout = (3, 1))
plot(tank_1_plot, tank_2_plot, layout = (2, 1))

# network_mn_copy = deepcopy(network_mn);
# WM._IM.update_data!(network_mn_copy, result_heur["solution"]);

# network_mn_copy = deepcopy(network_mn);
# WM._IM.update_data!(network_mn_copy, result_heur["solution"]);
# WM.set_start_all!(network_mn_copy);
# #WM.turn_on_all_components!(network_mn_copy);
# ext = Dict{Symbol, Any}(:pipe_breakpoints => 10, :pump_breakpoints => 10);
# result = WM.solve_mn_owf(network_mn_copy, WM.PWLRDWaterModel, gurobi_short; ext = ext);


# # Reload the network data.
# network = WM.parse_file(network_path);
# result = WaterModelsAnnex.simulate!(network, schedules, sched_sol, result_rel, ipopt);

# # Solve a relaxation of the OWF problem.
# ext = Dict{Symbol, Any}(:pipe_breakpoints => 10, :pump_breakpoints => 10);
# #result = WM.solve_mn_owf(network_mn, WM.PWLRDWaterModel, gurobi_2; ext = ext);
# result = solve_mn_owfh(network_mn, WM.PWLRDWaterModel, gurobi_2; ext = ext);

# ipopt = JuMP.optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0, "acceptable_tol" => 1.0e-8, "sb" => "yes");
# result = WaterModelsAnnex.simulate!(network, result, ipopt);
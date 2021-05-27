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

tank_1_max_vol = ones(24) * 0.25 * pi * network_heur["tank"]["1"]["diameter"]^2 * network_heur["tank"]["1"]["max_level"]
tank_1_milp = [result_milp["solution"]["nw"][string(n)]["tank"]["1"]["V"] for n in 1:24];
tank_1_milp_sim = [result_milp_sim["solution"]["nw"][string(n)]["tank"]["1"]["V"] for n in 1:24];
tank_1_fixed_sim = [result_heur_sim["solution"]["nw"][string(n)]["tank"]["1"]["V"] for n in 1:24];
tank_1_plot_milp = plot(1:6, [tank_1_milp[1:6], tank_1_milp_sim[1:6], tank_1_max_vol[1:6]], xlabel = "Time (h)", ylabel = "Tank 1 Volume", markershape = :circle, markersize = 3, label = ["MILP" "MILP Sim" "Max Volume"], legend=(0.2, 0.4));
tank_1_plot = plot(1:24, [tank_1_milp, tank_1_fixed_sim], xlabel = "Time (h)", ylabel = "Tank 1 Volume", markershape = :circle, markersize = 3, label = ["MILP" "Repaired"], legendfontsize = 8, legend=(0.15, 0.35));

tank_2_max_vol = ones(24) * 0.25 * pi * network_heur["tank"]["2"]["diameter"]^2 * network_heur["tank"]["2"]["max_level"]
tank_2_milp = [result_milp["solution"]["nw"][string(n)]["tank"]["2"]["V"] for n in 1:24];
tank_2_milp_sim = [result_milp_sim["solution"]["nw"][string(n)]["tank"]["2"]["V"] for n in 1:24];
tank_2_fixed_sim = [result_heur_sim["solution"]["nw"][string(n)]["tank"]["2"]["V"] for n in 1:24];
tank_2_plot_milp = plot(1:6, [tank_2_milp[1:6], tank_2_milp_sim[1:6], tank_2_max_vol[1:6]], xlabel = "Time (h)", ylabel = "Tank 2 Volume", markershape = :circle, markersize = 3, legend = false);
tank_2_plot = plot(1:24, [tank_2_milp, tank_2_fixed_sim], xlabel = "Time (h)", ylabel = "Tank 2 Volume", markershape = :circle, markersize = 3, legend = false);

pump_1_milp = [result_milp["solution"]["nw"][string(n)]["pump"]["1"]["status"] for n in 1:24];
pump_1_fixed_sim = [Int(result_heur_sim["solution"]["nw"][string(n)]["pump"]["1"]["status"]) for n in 1:24];
pump_1_plot = plot(1:24, [pump_1_milp, pump_1_fixed_sim], xlabel = "Time (h)", ylabel = "Pump 1 Status", markershape = :circle, markersize = 3, label = ["MILP" "Repaired"], legendfontsize = 8, legend = (0.83, 0.5));

pump_2_milp = [result_milp["solution"]["nw"][string(n)]["pump"]["2"]["status"] for n in 1:24];
pump_2_fixed_sim = [Int(result_heur_sim["solution"]["nw"][string(n)]["pump"]["2"]["status"]) for n in 1:24];
pump_2_plot = plot(1:24, [pump_2_milp, pump_2_fixed_sim], xlabel = "Time (h)", ylabel = "Pump 2 Status", markershape = :circle, markersize = 3, legend = false);

pump_3_milp = [result_milp["solution"]["nw"][string(n)]["pump"]["3"]["status"] for n in 1:24];
pump_3_fixed_sim = [Int(result_heur_sim["solution"]["nw"][string(n)]["pump"]["3"]["status"]) for n in 1:24];
pump_3_plot = plot(1:24, [pump_3_milp, pump_3_fixed_sim], xlabel = "Time (h)", ylabel = "Pump 3 Status", markershape = :circle, markersize = 3, legend = false);


plot(tank_1_plot_milp, tank_2_plot_milp, layout = (2, 1))
plot(tank_1_plot, tank_2_plot, layout = (2, 1))
plot(pump_1_plot, pump_2_plot, pump_3_plot, layout = (3, 1), legendfontsize = 4)
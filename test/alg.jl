@testset "Network Expansion Algorithms" begin
    @testset "Shamir network, No-good Cutting Plane MILP-R Method." begin
        network_path = "../test/data/epanet/shamir.inp"
        modifications_path = "../test/data/json/shamir.json"
        cut_method = WaterModelsAnnex.add_no_good_cut!
        WaterModelsAnnex.ne_global(network_path, modifications_path, ipopt, cbc, cut_method, "")
    end
end

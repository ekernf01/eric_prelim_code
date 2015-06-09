workspace()
push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")
using Chem_rxn_tools
using Distributions
using Winston

t_obs = [1:50.0]
SBML_file = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/wilkinson_rxns_SBML_shorthand.txt"
wilk_cri = Chem_rxn_tools.SBML_read(SBML_file)
obs_mol_name = "SigD"
mols_to_show = ["SigD", "CodY"]
noise_distribution = Normal(0,10)
sim_results1 = Chem_rxn_tools.make_sim_data(t_obs, wilk_cri, obs_mol_name, noise_distribution)
sim_results2 = Chem_rxn_tools.gillespie(wilk_cri, t_obs[end], false)
test_plot_folder = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/Chem_rxn_tools/test"
Chem_rxn_tools.plot_save_sim_data(test_plot_folder, sim_results1, wilk_cri, mols_to_show)
n = sim_results2.num_rxns_occ
obs_vecs = [[sim_results2.x_path[i][j]::Int64 for i in 1:sim_results2.num_rxns_occ] for j in 1:length(sim_results2.x_path[1])]
plot(sim_results2.rxn_times, obs_vecs[2]-obs_vecs[1])
wilk_cri.sto_mat_nonzero_inds

println(wilk_cri.init_amts)

obs_vecs[4]
println(sim_results2.x_path[10])
println(sim_results2.x_path[10])
println(sim_results2.x_path[10])


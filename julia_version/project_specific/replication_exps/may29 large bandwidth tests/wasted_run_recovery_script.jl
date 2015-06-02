workspace()
push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")
using Chem_rxn_tools
  SBML_file = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/wilkinson_rxns_SBML_shorthand.txt"
  wilk_cri = Chem_rxn_tools.SBML_read(SBML_file)
  sto_mat_graph,rxn_entry_mat_graph = Chem_rxn_tools.make_cri_graphic(wilk_cri)
using pMCMC_julia
using HDF5, JLD
cd()
today_filepath = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/replication_exps/may29 large bandwidth tests/wasted_run"
cd(today_filepath)
readdir()

MCS = pMCMC_julia.MCMC_state()

@load "stage23sample"
MCS.current_sample = current_sample
unk_names = ["SigDprod", "flacherep", "flacheunrep"]
unk_rates = [1.0, 0.02, 0.1]
unk_inds = Int64[]
for word in unk_names
  push!(unk_inds, Chem_rxn_tools.get_rate_indices(wilk_cri, word))
end
include("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/contour_bivariate_plot_maker.jl")
post_contour = contour_plot_two_mols(MCS, today_filepath, unk_names[2], unk_names[3], wilk_cri, unk_rates, unk_names)
include("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/bivariate_plot_maker.jl")
post_biv = plot_save_two_mols(MCS, today_filepath, unk_names[2], unk_names[3], wilk_cri, unk_rates, unk_names)


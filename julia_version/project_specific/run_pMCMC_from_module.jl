#This allows Julia to find Eric's modules.
push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")

#Load in info from SBML shorthand and visualize results.
using Chem_rxn_tools
  SBML_file = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/wilkinson_rxns_SBML_shorthand.txt"
  wilk_cri = Chem_rxn_tools.SBML_read(SBML_file)
  sto_mat_graph,rxn_entry_mat_graph = Chem_rxn_tools.make_cri_graphic(wilk_cri)
sto_mat_graph
rxn_entry_mat_graph


using pMCMC_julia
  prior_sample =
  current_MCMC_state = MCMC_state()

pMCMC(d_obs, t_obs, current_MCMC_state, emission_logden, fwd_sim)




workspace()
  #This allows Julia to find Eric's modules.
  push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")
  #Load in info from SBML shorthand and visualize results.
  using Chem_rxn_tools
  SBML_file = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/wilkinson_rxns_SBML_shorthand.txt"
  wilk_cri = Chem_rxn_tools.SBML_read(SBML_file)
  sto_mat_graph,rxn_entry_mat_graph = Chem_rxn_tools.make_cri_graphic(wilk_cri)

Chem_rxn_tools.gillespie_basic_test()
decay_plot = Chem_rxn_tools.gillespie_stoch_mean_test(20, 0.01*ones(Float64, 1), 100)
decay_plot = Chem_rxn_tools.gillespie_stoch_mean_test_sampler_version(200, 0.01*ones(Float64, 1), 100)

sto_mat_graph
rxn_entry_mat_graph


#This saves the plots.
#using Winston
#  savefig(sto_mat_graph, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/sto_mat_graph.png","width",1536,"height",360)
#  savefig(rxn_entry_mat_graph, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/rxn_entry_mat_graph.png","width",1536,"height",360)

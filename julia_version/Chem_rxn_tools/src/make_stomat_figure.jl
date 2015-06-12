#This allows Julia to find Eric's modules.
workspace()
  push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")

#-----------------------------Load in info from SBML shorthand and visualize results.------------------------------
  using Chem_rxn_tools
  using Winston
  SBML_file = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/wilkinson_rxns_SBML_shorthand.txt"
  wilk_cri = Chem_rxn_tools.SBML_read(SBML_file)
  sto_mat_graph,rxn_entry_mat_graph = Chem_rxn_tools.make_cri_graphic(wilk_cri)
  savefig(sto_mat_graph,"/Users/EricKernfeld/Desktop/sto_fig.png","width",1500,"height",450)
  savefig(rxn_entry_mat_graph,"/Users/EricKernfeld/Desktop/rmat_fig_.png","width",1500,"height",450)



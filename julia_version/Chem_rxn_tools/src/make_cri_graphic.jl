function make_cri_graphic(cri::Chem_rxn_info)
  rxn_rate_graph_labels = String[]
  for rate in cri.rxn_rates
    push!(rxn_rate_graph_labels, string(rate))
  end

  #-------------------------------Plot matrices with stoichiometric details-------------------------------
  sto_mat_graph = imagesc(cri.sto_mat')
    setattr(sto_mat_graph.x1,
            label="Molecules",
            ticks=(1:cri.num_species)-0.5,
            draw_ticks=false,
            ticklabels=cri.species_labels)
    setattr(sto_mat_graph.y1,
            label="Stoichiometry",
            ticks=(1:cri.num_rxns)-0.5,
            draw_ticks=false,
            ticklabels=cri.rxns_written_out)
    sto_mat_graph
    setattr(sto_mat_graph.y2,
            label="Rates",
            ticks=(1:cri.num_rxns)-0.5,
            draw_ticks=false,
            ticklabels=rxn_rate_graph_labels)

  rxn_entry_mat_graph = imagesc(cri.rxn_entry_mat')
    setattr(rxn_entry_mat_graph.x1,
            label="Molecules",
            ticks=(1:cri.num_species)-0.5,
            draw_ticks=false,
            ticklabels=cri.species_labels)
    setattr(rxn_entry_mat_graph.y1,
            label="Stoichiometry",
            ticks=(1:cri.num_rxns)-0.5,
            draw_ticks=false,
            ticklabels=cri.rxns_written_out)
    setattr(rxn_entry_mat_graph.y2,
            label="Rates",
            ticks=(1:cri.num_rxns)-0.5,
            draw_ticks=false,
            ticklabels=rxn_rate_graph_labels)

  return sto_mat_graph , rxn_entry_mat_graph
end

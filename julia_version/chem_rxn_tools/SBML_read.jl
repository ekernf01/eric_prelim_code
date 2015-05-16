
#Reads SBML shorthand, turning it into input for gillespie().
#This is a primitive module. It currently ignores everything but species, parameters and reactions.
# It will also make errors if one species name contains another, like in this reaction.
# CodYPhaggfp->CodY+Phag_gfp
#
#I used a cheap work-around that looks for _, so that CodYPhag_gfp will be read as
# if CodY is present, but CodY_Phag_gfp will not have the problem. This works well with
# Darren Wilkinson's naming conventions.
#
#There are probably other undetected problems. There is no guarantee that this script can read SBML shorthand even close to properly.
#
#output:
# init_x, the particle counts
# sto_mat, the stoichiometry matrix, the matrix whose i,j entry says how many
#      molecules of type i are consumed by a rxn of type j (net change)
# rxn_entry_mat, the matrix whose i,j entry says how many molecules
#      of type i enter a rxn of type j (not a net change)
# rxn_rates--vector of floats. rates of reactions encoded, MEASURED IN INTENSITY PER SECOND.
# num_species--scalar integer. how many different molecules are encoded
# species_labels--vector of strings. names of the molecules encoded
# This list of output is out of date

function SBML_read(SBML_file)

  wilk_SBML_file = open(SBML_file)
  lines = readlines(wilk_SBML_file)

  wilk_chem_rxn_info = chem_rxn_info()

  get_par_info!(lines,wilk_chem_rxn_info)
  get_species_info!(lines,wilk_chem_rxn_info)
  get_rxn_info!(lines,wilk_chem_rxn_info)

  println("GFP is at index ", chem_indices.GFP_ind,
          " and species_info.species_labels[", chem_indices.GFP_ind, "] is ", species_labels[chem_indices.GFP_ind], ".")
  println("SigD is at index ", chem_indices.SigD_ind,
          " and species_info.species_labels[", chem_indices.SigD_ind, "] is ", species_labels[chem_indices.SigD_ind], ".")
  println("Hag is at index ", chem_indices.Hag_ind,
          " and species_info.species_labels[", chem_indices.Hag_ind, "] is ", species_labels[chem_indices.Hag_ind], ".")
  println("CodY is at index ", chem_indices.CodY_ind,
          " and species_info.species_labels[", chem_indices.CodY_ind, "] is ", species_labels[chem_indices.CodY_ind], ".")

  #-------------------------------Fill in sto_mat and rxn_entry_mat-------------------------------


#-------------------------------  #Fill in true_rxn_rates-------------------------------
  num_params_at_start = length(par_names)
  true_rxn_rates = zeros(Float64,num_rxns)
  for j in 1:length(rxn_pos_in_file)
    rxn_pos = rxn_pos_in_file[j]
    line_with_param = strip(lines[rxn_pos + 2])
    for i in 1:num_params_at_start
      #if there's an actual number, it's after the =. Grab it.
      if contains(line_with_param,"=")
        eq_ind = findin(line_with_param, "=")[1]
        true_rxn_rates[j] = float(line_with_param[(eq_ind + 1):end])
      end
      #if there's no actual number, find the parameter that was set earlier.
      if contains(line_with_param,par_names[i])
        true_rxn_rates[j] = par_vals[i]
      end
    end
  end
  true_rxn_rates

  rxn_rate_graph_labels = String[]
  for rate in true_rxn_rates
    push!(rxn_rate_graph_labels, string(rate))
  end
  rxn_rate_graph_labels



#-------------------------------Plot matrices with stoichiometric details-------------------------------
  sto_mat_graph = imagesc(sto_mat')
    setattr(sto_mat_graph.x1,
            label="Molecules",
            ticks=(1:num_species)-0.5,
            draw_ticks=false,
            ticklabels=species_labels)
    setattr(sto_mat_graph.y1,
            label="Stoichiometry",
            ticks=(1:num_rxns)-0.5,
            draw_ticks=false,
            ticklabels=rxns_written_out)
    setattr(sto_mat_graph.y2,
            label="Rates",
            ticks=(1:num_rxns)-0.5,
            draw_ticks=false,
            ticklabels=rxn_rate_graph_labels)
    sto_mat_graph
    savefig(sto_mat_graph, "sto_mat_graph.png","width",1536,"height",360)

  rxn_entry_mat_graph = imagesc(rxn_entry_mat')
    setattr(rxn_entry_mat_graph.x1,
            label="Molecules",
            ticks=(1:num_species)-0.5,
            draw_ticks=false,
            ticklabels=species_labels)
    setattr(rxn_entry_mat_graph.y1,
            label="Stoichiometry",
            ticks=(1:num_rxns)-0.5,
            draw_ticks=false,
            ticklabels=rxns_written_out)
    setattr(rxn_entry_mat_graph.y2,
            label="Rates",
            ticks=(1:num_rxns)-0.5,
            draw_ticks=false,
            ticklabels=rxn_rate_graph_labels)
    rxn_entry_mat_graph
    savefig(rxn_entry_mat_graph, "rxn_entry_mat_graph.png","width",1536,"height",360)
  return init_x, sto_mat, rxn_entry_mat, true_rxn_rates, num_species, species_labels, num_rxns, rxn_labels, GFP_ind, SigD_ind, Hag_ind, CodY_ind
end


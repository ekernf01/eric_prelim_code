
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

function SBML_read(SBML_file_path)
  println("Warning: SBML_read is not guaranteed to work properly. Check your results. Maybe make_cri_graphic() would help.")
  SBML_file = open(SBML_file_path)
  lines = readlines(SBML_file)

  my_cri = Chem_rxn_info()

  get_par_info!(lines,my_cri)
  get_species_info!(lines,my_cri)
  get_rxn_info!(lines,my_cri)
  get_stoich_info!(lines,my_cri)
  get_rate_info!(lines,my_cri)
  chem_rxn_data_check!(my_cri)
#   println("GFP is at index ", chem_indices.GFP_ind,
#           " and species_info.species_labels[", chem_indices.GFP_ind, "] is ", species_labels[chem_indices.GFP_ind], ".")
#   println("SigD is at index ", chem_indices.SigD_ind,
#           " and species_info.species_labels[", chem_indices.SigD_ind, "] is ", species_labels[chem_indices.SigD_ind], ".")
#   println("Hag is at index ", chem_indices.Hag_ind,
#           " and species_info.species_labels[", chem_indices.Hag_ind, "] is ", species_labels[chem_indices.Hag_ind], ".")
#   println("CodY is at index ", chem_indices.CodY_ind,
#           " and species_info.species_labels[", chem_indices.CodY_ind, "] is ", species_labels[chem_indices.CodY_ind], ".")

  return my_cri
end


#-------------------------------Fill in sto_mat and rxn_entry_mat-------------------------------
function get_stoich_info!(lines, cri::Chem_rxn_info)
  ns = cri.num_species
  nr = cri.num_rxns
  cri.sto_mat = zeros(Int64, ns, nr)
  cri.rxn_entry_mat = zeros(Int64, ns, nr)

  for j in 1:nr
    rxn_pos = cri.rxn_pos_in_SBML_file[j]
    line_with_stoich = strip(lines[rxn_pos + 1])
    arrow_ind = findin(line_with_stoich, "->")
    rxn_LHS = line_with_stoich[1:(arrow_ind[1]-1)]
    rxn_RHS = line_with_stoich[(arrow_ind[2]+1):end]
    for i in 1:ns
      present_on_LHS = contains(rxn_LHS, cri.species_labels[i])&
          !contains(rxn_LHS, string(cri.species_labels[i],"_"))&
          !contains(rxn_LHS, string("_",cri.species_labels[i]))
      present_on_RHS = contains(rxn_RHS, cri.species_labels[i])&
          !contains(rxn_RHS, string(cri.species_labels[i],"_"))&
          !contains(rxn_RHS, string("_",cri.species_labels[i]))
      if present_on_LHS
        cri.sto_mat[i,j] = cri.sto_mat[i,j] - 1
        cri.rxn_entry_mat[i,j] = cri.rxn_entry_mat[i,j] + 1
        #println(string(rxn_LHS, " contains ", species_labels[i], " so decrement sto_mat[", i,",", j,"] by one."))
      end
      if present_on_RHS
        cri.sto_mat[i,j] = cri.sto_mat[i,j] + 1
        #println(string(rxn_RHS, " contains ", species_labels[i], " so increment sto_mat[", i,",", j,"] by one."))
      end
    end
    push!(cri.rxns_written_out, line_with_stoich)
  end
end


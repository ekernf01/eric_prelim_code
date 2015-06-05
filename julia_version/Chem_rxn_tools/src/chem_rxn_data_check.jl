
# ——————————————————————————————————————————————
function chem_rxn_data_check!(cri_to_check::Chem_rxn_info)
# rxn_entry_mat is same size as sto_mat
# length of init_x matches num rows of sto_mat
# length of rxn_rates matches num cols of sto_mat
# init_x, sto_mat, rxn_entry_mat are nonnegative integers

  if size(cri_to_check.rxn_entry_mat)!=size(cri_to_check.sto_mat)
    error("chem_rxn_info object has matrices unmatched in size.")
    return false
  end
  if size(cri_to_check.sto_mat,1)!=size(cri_to_check.init_amts,1)
    error("chem_rxn_info object has sto_mat & init_x unmatched in size.")
    return false
  end
  if size(cri_to_check.sto_mat, 2)!=size(cri_to_check.rxn_rates,1)
    println(size(cri_to_check.sto_mat))
    println(size(cri_to_check.rxn_rates))
    error("chem_rxn_info object has sto_mat & rxn_rates unmatched in size.")
    return false
  end
  if sum(cri_to_check.init_amts.<0)>0
    error("chem_rxn_info object has negative init_x")
    return false
  end
  #Negative values should be ok here!
#   if sum(sto_mat.<0)>0
#     error("chem_rxn_info object has negative sto_mat")
#     return false
#   end
  if sum(cri_to_check.rxn_entry_mat .< 0)>0
    error("chem_rxn_info object has negative rxn_entry_mat")
    return false
  end

  cri_to_check.sto_mat_nonzero_inds = Array{Int64,1}[]
  cri_to_check.rxn_entry_mat_nonzero_inds = Array{Int64,1}[]
  for rxn_ind in 1:cri_to_check.num_rxns
    rxn_entry_temp_nonzeros = Int64[];
    sto_mat_temp_nonzeros = Int64[];
    for mol_ind in 1:cri_to_check.num_species
      if cri_to_check.sto_mat[mol_ind,rxn_ind]!=0
        push!(sto_mat_temp_nonzeros, mol_ind)
      end
      if cri_to_check.rxn_entry_mat[mol_ind,rxn_ind]!=0
        push!(rxn_entry_temp_nonzeros, mol_ind)
      end
    end
    push!(cri_to_check.sto_mat_nonzero_inds, sto_mat_temp_nonzeros)
    push!(cri_to_check.rxn_entry_mat_nonzero_inds, rxn_entry_temp_nonzeros)
  end
  println("chem_rxn_data_check! returning with cri_to_check.sto_mat_nonzero_inds = ", cri_to_check.sto_mat_nonzero_inds)
   println("chem_rxn_data_check! returning with cri_to_check.rxn_entry_mat_nonzero_inds = ", cri_to_check.rxn_entry_mat_nonzero_inds)
  return true
end


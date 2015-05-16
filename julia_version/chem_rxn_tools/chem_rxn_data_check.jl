
# ——————————————————————————————————————————————
function chem_rxn_data_check(cri_to_check::Chem_rxn_info)
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
  return true
end


#---------------------------------Automatically finds which index holds the reaction rates we'll treat as unknown---------------------------------
#Give it a string and it finds the chemical with that name
function get_rate_indices(cri::Chem_rxn_info, rxn_name::String)
  for i in 1:cri.num_rxns
    word = cri.rxn_labels[i]
    rxn_present = contains(word, rxn_name)&
          !contains(word, string(rxn_name,"_"))&
          !contains(word, string("_",rxn_name))
    if rxn_present
      return  i
    end
  end
  return 0
end

#Unit test
function get_rate_indices_tester()
  demo_cri = Chem_rxn_info()
  demo_cri.rxn_labels = ["kryptonite", "kryptonite_cinnamon_rolls", "cinnamon_rolls"]
  demo_cri.num_rxns = 3
  chem_rxn_data_check(demo_cri)
  flag1 = (1==get_rate_indices(demo_cri,"kryptonite"))

  demo_cri.rxn_labels = ["kryptonite_cinnamon_rolls", "cinnamon_rolls"]
  demo_cri.num_rxns = 2
  chem_rxn_data_check(demo_cri)
  flag2 = (0==get_rate_indices(demo_cri, "kryptonite"))

  if flag1 && flag2
    return 0
  else
    println(flag1)
    println(flag2)
    error("get_chem_indices is not correct.")
  end
end

get_rate_indices_tester()



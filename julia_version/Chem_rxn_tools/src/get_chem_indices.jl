
#---------------------------------Automatically finds which index holds the chemicals we'll measure or plot---------------------------------
#Give it a string and it finds the chemical with that name
function get_chem_indices(cri::Chem_rxn_info, mol_name::String)
  for i in 1:cri.num_species
    word = cri.species_labels[i]
    mol_present = contains(word, mol_name)&
          !contains(word, string(mol_name,"_"))&
          !contains(word, string("_",mol_name))
    if mol_present
      return  i
    end
  end
  return 0
end

#Unit test
function get_chem_indices_tester()
  demo_cri = Chem_rxn_info()
  demo_cri.species_labels = ["kryptonite", "kryptonite_cinnamon_rolls", "cinnamon_rolls"]
  demo_cri.num_species = 3
  chem_rxn_data_check!(demo_cri)
  flag1 = (1==get_chem_indices(demo_cri,"kryptonite"))

  demo_cri.species_labels = ["kryptonite_cinnamon_rolls", "cinnamon_rolls"]
  demo_cri.num_species = 2
  chem_rxn_data_check!(demo_cri)
  flag2 = (0==get_chem_indices(demo_cri, "kryptonite"))

  if flag1 && flag2
    return 0
  else
    println(flag1)
    println(flag2)
    error("get_chem_indices is not correct.")
  end
end

get_chem_indices_tester()



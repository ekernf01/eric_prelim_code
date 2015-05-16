

#---------------------------------fill in positions of reaction in SBML file, reaction names, and how many--------------------------------
function get_rxn_info!(lines::Array{String, 1}, cri::Chem_rxn_info)
  for i in 1:length(lines)
    if contains(lines[i],"@r=") #If you see SBML for a new reaction,
      cri.num_rxns = cri.num_rxns + 1
      push!(cri.rxn_labels, replace(lines[i], "@r=", ""))
      push!(cri.rxn_pos_in_SBML_file,i)
    end
  end
end

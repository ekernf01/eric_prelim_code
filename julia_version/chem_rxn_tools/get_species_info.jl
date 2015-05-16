#---------------------------------fill in species_labels, init_x and num_species---------------------------------
#Parse the molecule types.
function get_species_info!(lines::Array{String, 1}, cri::Chem_rxn_info)
  cri.num_species = 1
  for i in 1:length(lines)
    #If you see the section where it defines the species
    #Loop over lines until you hit a new SBML code block
    if contains(lines[i],"@species")
      while !contains(lines[i+cri.num_species],"@")
        colon_in = findin(lines[i+cri.num_species], ":")[1]
        eq_ind = findin(lines[i+cri.num_species], "=")[1]
        push!(cri.species_labels, lines[i+cri.num_species][(colon_in+1):(eq_ind-1)]) #these parens are actually redundant in Julia.
        init_x_str = lines[i+num_species][(eq_ind+1):end]

        #remove some extraneous s's. Not sure what they do in original SBML.
        push!(cri.init_amts, float(replace(init_x_str, "s", "")))
        cri.num_species = cri.num_species + 1
      end
      species_info.num_species = species_info.num_species - 1 #corrects an off-by-one
    end
  end
end

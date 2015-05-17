

#---------------------------------fill in rxn rate names and values---------------------------------
function get_par_info!(lines, cri::Chem_rxn_info)
  for i in 1:length(lines)
    if contains(lines[i],"@parameters")
      j = i + 1
      while(!contains(lines[j], "@")) #Loop over lines until you hit a new SBML code block
        eq_ind = findin(lines[j], "=")[1]
        push!(cri.SBML_par_names, strip(lines[j][1:(eq_ind-1)])) #remove whitespace
        push!(cri.SBML_par_vals, float(lines[j][(eq_ind+1):end])) #float discards whitespace automatically
        j = j + 1
      end
    end
  end
end


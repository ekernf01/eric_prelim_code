#Reads SBML shorthand, turning it into input for gillespie().
#This is a primitive script. It currently ignores everything but species, parameters and reactions.
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
function SBML_read()

  cd(string("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/Transcribing_SBML/"))
  wilk_SBML_file = open("wilkinson_rxns_SBML_shorthand.txt")
  lines = readlines(wilk_SBML_file)

  #---------------------------------fill in par_names and par_vals---------------------------------
  par_names =  String[]
  par_vals = Float64[]
  for i in 1:length(lines)
    if contains(lines[i],"@parameters")
      j = i + 1
      while(!contains(lines[j], "@")) #Loop over lines until you hit a new SBML code block
        eq_ind = findin(lines[j], "=")[1]
        push!(par_names, strip(lines[j][1:(eq_ind-1)])) #remove whitespace
        push!(par_vals, float(lines[j][(eq_ind+1):end])) #float discards whitespace automatically
        j = j + 1
      end
    end
  end

  #---------------------------------fill in species_labels, init_x and num_species---------------------------------
  #Parse the molecule types.
  species_labels = String[]
  init_x = Int64[]
  num_species = 1
  for i in 1:length(lines)
    if contains(lines[i],"@species") #If you see the section where it defines the species
      while !contains(lines[i+num_species],"@")#Loop over lines until you hit a new SBML code block
        colon_in = findin(lines[i+num_species], ":")[1]
        eq_ind = findin(lines[i+num_species], "=")[1]
        push!(species_labels, lines[i+num_species][(colon_in+1):(eq_ind-1)]) #these parens are actually redundant in Julia.
        init_x_str = lines[i+num_species][(eq_ind+1):end]
        push!(init_x, float(replace(init_x_str, "s", ""))) #remove some extraneous s's. Not sure what this does in original SBML.
        num_species = num_species + 1
      end
      num_species = num_species - 1 #corrects an off-by-one
    end
  end
  species_labels
  init_x
  num_species


  #---------------------------------fill in rxn_pos_in_file, rxn_labels, num_rxns---------------------------------
  rxn_labels = String[]
  num_rxns = 0
  rxn_pos_in_file = Int64[]
  for i in 1:length(lines)
    if contains(lines[i],"@r=") #If you see SBML for a new reaction,
      num_rxns = num_rxns + 1
      push!(rxn_labels, replace(lines[i], "@r=", ""))
      push!(rxn_pos_in_file,i)
    end
  end
  rxn_labels
  num_rxns
  rxn_pos_in_file

  #-------------------------------Fill in sto_mat and rxn_entry_mat-------------------------------
  sto_mat = zeros(num_species, num_rxns)
  rxn_entry_mat = zeros(num_species, num_rxns)
  rxns_written_out = String[]
  for j in 1:length(rxn_pos_in_file)
    rxn_pos = rxn_pos_in_file[j]
    line_with_stoich = strip(lines[rxn_pos + 1])
    arrow_ind = findin(line_with_stoich, "->")
    rxn_LHS = line_with_stoich[1:(arrow_ind[1]-1)]
    rxn_RHS = line_with_stoich[(arrow_ind[2]+1):end]
    for i in 1:num_species
      present_on_LHS = contains(rxn_LHS, species_labels[i])&
          !contains(rxn_LHS, string(species_labels[i],"_"))&
          !contains(rxn_LHS, string("_",species_labels[i]))
      present_on_RHS = contains(rxn_RHS, species_labels[i])&
          !contains(rxn_RHS, string(species_labels[i],"_"))&
          !contains(rxn_RHS, string("_",species_labels[i]))
      if present_on_LHS
        sto_mat[i,j] = sto_mat[i,j] - 1
        rxn_entry_mat[i,j] = rxn_entry_mat[i,j] + 1
        #println(string(rxn_LHS, " contains ", species_labels[i], " so decrement sto_mat[", i,",", j,"] by one."))
      end
      if present_on_RHS
        sto_mat[i,j] = sto_mat[i,j] + 1
        #println(string(rxn_RHS, " contains ", species_labels[i], " so increment sto_mat[", i,",", j,"] by one."))
      end
    end
    push!(rxns_written_out, line_with_stoich)
  end



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

  sto_mat
  rxn_entry_mat
  rxns_written_out
  using Winston

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

  return init_x, sto_mat, rxn_entry_mat, true_rxn_rates, num_species, species_labels, num_rxns, rxn_labels, GFP_ind, sigD_ind
end


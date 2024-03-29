#Given a Chem_rxn_info object and a simulation time (a positive real scalar), gillespie returns either:
#  if inside_sampler is true, the molecule counts at time T_sim
#  if inside_sampler is false, a Chem_sim_result object.

function gillespie(cri::Chem_rxn_info, T_sim::Float64, inside_sampler::Bool)

  current_x = copy(cri.init_amts)
  num_rxns_occ = 0
  t_spent = 0

  if(!inside_sampler)
    if(T_sim<0 || !chem_rxn_data_check!(cri))
      error("gillespie received bad input.")
    end
    rxn_times = Array(Float64, 0)
    rxn_types = Array(Int64, 0)
    x_path = Array{Int64,1}[]
    push!(x_path, cri.init_amts)
    push!(rxn_times, 0.0)
  end
  unit_rate_expo = Exponential(1)
  while(t_spent < T_sim)
    #Get reaction propensities, prod_j c_i* (X_j choose k_ij)
    alpha = copy(cri.rxn_rates)

    # For the jth reaction, sto_mat_nonzero_inds[j] is an array of indices so
    #that sto_mat[sto_mat_nonzero_inds[rxn_index][i], rxn_index] is nonzero (and nothing else is).
    if !isempty(cri.rxn_entry_mat_nonzero_inds)
      for rxn_index = 1:cri.num_rxns
        for mol_index in cri.rxn_entry_mat_nonzero_inds[rxn_index]
          alpha[rxn_index] = alpha[rxn_index]*binomial(current_x[mol_index], cri.rxn_entry_mat[mol_index, rxn_index])
        end
      end
    end
    alpha_sum = sum(alpha)

    #Perchance we've no propensity to react, stop the simulation.
    #(It'll stop: see next block)
    #Otherwise, step forward in time
    if alpha_sum == 0
      t_spent = T_sim
    else
      tau = rand(unit_rate_expo)/alpha_sum
      t_spent = t_spent + tau
    end

    #If time's not up, update the molecule counts, the reactions
    #count, and the reaction times and types list
    if(t_spent < T_sim)
      current_rxn_type = rand(Categorical(alpha/alpha_sum))
      if !isempty(cri.sto_mat_nonzero_inds)
        for mol_index in cri.sto_mat_nonzero_inds[current_rxn_type]
          current_x[mol_index] = current_x[mol_index] + cri.sto_mat[mol_index,current_rxn_type]
        end
      end

      #! means that this modifies the array that it is given
      #also means logical negation
      if(!inside_sampler)
        #println("In gillespie, current_x: ", current_x, " t: ", t_spent, " type: ", current_rxn_type)
        num_rxns_occ = num_rxns_occ + 1
        push!(rxn_times, t_spent)
        push!(rxn_types, current_rxn_type)
        push!(x_path, copy(current_x))
      end
    end
  end
  if(inside_sampler)
    return current_x
  else
    x_obs = Array{Int64, 1}[]
    d_obs = Float64[]
    t_obs = Float64[]
    obs_mol_name = ""
    obs_mol_ind = -1

    return Chem_sim_result(x_path,
                           current_x,
                           num_rxns_occ,
                           rxn_types,
                           rxn_times,
                           t_spent,
                           x_obs,
                           d_obs,
                           t_obs,
                           obs_mol_name,
                           obs_mol_ind)
  end
end

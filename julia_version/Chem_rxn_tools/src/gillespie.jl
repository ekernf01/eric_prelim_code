
# ——————————————————————————————————————————————
# Forward simulation files:
# gillespie.jl

# At the core, written partly in C and called by many other parts of the project, is the forward simulator.

# Unit tests should:
# --feed it garbage: strings, negative numbers and fractions where only positive numbers or integers should be, matrices where it needs vectors, stuff that's the wrong length
# --Run it many times on a simple reaction model and check to see it converge to the true stochastic mean
# ——————————————————————————————————————————————
# Given:

# init_x, the particle counts
# sto_mat, the stoichiometry matrix, the matrix whose i,j entry says how many
#      molecules of type i are consumed by a rxn of type j (net change)
# rxn_entry_mat, the matrix whose i,j entry says how many molecules
#      of type i enter a rxn of type j (not a net change)
# rxn_rates, reaction rates, MEASURED IN INTENSITY PER SECOND
# T_sim, the time over which to simulate, MEASURED IN SECONDS
# inside_sampler, a boolean telling it whether this run is inside the LF-pMCMC
# sto_mat_nonzero_inds and rxn_entry_mat_nonzero_inds are arrays that help take advantage of sparsity.
# For the jth reaction, sto_mat_nonzero_inds[j] is an array of indices so that sto_mat[sto_mat_nonzero_inds[j][i], j] is nonzero (and nothing else is).
# Similar for rxn_entry_mat.

# It verifies:

# rxn_entry_mat is same size as sto_mat
# init_x, sto_mat, rxn_entry_mat are nonnegative integers
# length of init_x matches num rows of sto_mat
# length of rxn_rates matches num cols of sto_mat
# T_sim is a nonnegative float or double

# Then, it runs the Gillespie algorithm to produce:

# rxn_times, vector of times at which reactions occur (for convenience, rxn_time[1] should be 0)
# rxn_types, list of reactions that occured (integers that index sto_mat)
# x_path, a list of vectors with first element init_x and ith element showing the molecule counts between rxn_time[i-1] and rxn_time[i].

function gillespie(cri::Chem_rxn_info, T_sim::Float64, inside_sampler::Bool)

  current_x = cri.init_amts
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

    # For the jth reaction, sto_mat_nonzero_inds[j] is an array of indices so
    #that sto_mat[sto_mat_nonzero_inds[rxn_index][i], rxn_index] is nonzero (and nothing else is).
    alpha = copy(cri.rxn_rates)
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
          current_x = current_x + cri.sto_mat[mol_index,current_rxn_type]
        end
      end

      #! means that this modifies the array that it is given
      #also means logical negation
      if(!inside_sampler)
        num_rxns_occ = num_rxns_occ + 1
        push!(rxn_times, t_spent)
        push!(rxn_types, current_rxn_type)
        push!(x_path, current_x)
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


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

# It verifies:

# rxn_entry_mat is same size as sto_mat
# init_x, sto_mat, rxn_entry_mat are nonnegative integers
# length of init_x matches num rows of sto_mat
# length of rxn_rates matches num cols of sto_mat
# T_sim is a nonnegative float or double

# Then, it runs the Gillespie algorithm to produce:

# rxn_times, vector of times at which reactions occur (for convenience, rxn_time[0] should be 0)
# rxn_types, list of reactions that occured (integers that index sto_mat)
# x_path, a list of vectors with zero element init_x and i element showing the molecule counts between rxn_time[i] and rxn_time[i].

# ——————————————————————————————————————————————

function gillespie_input_ok(init_x, sto_mat,rxn_entry_mat, rxn_rates, T_sim)
# rxn_entry_mat is same size as sto_mat
# length of init_x matches num rows of sto_mat
# length of rxn_rates matches num cols of sto_mat
# init_x, sto_mat, rxn_entry_mat are nonnegative integers
# T_sim is a nonnegative float or double
  if size(rxn_entry_mat)!=size(sto_mat)
    error("gillespie received matrices unmatched in size.")
    return false
  end
  if size(sto_mat,1)!=size(init_x,1)
    error("gillespie received sto_mat & init_x unmatched in size.")
    return false
  end
  if size(sto_mat, 2)!=size(rxn_rates,1)
    error("gillespie received sto_mat & rxn_rates unmatched in size.")
    return false
  end
  if sum(init_x.<0)>0
    error("gillespie received negative init_x")
    return false
  end
  #Negative values should be ok here!
#   if sum(sto_mat.<0)>0
#     error("gillespie received negative sto_mat")
#     return false
#   end
  if sum(rxn_entry_mat.<0)>0
    error("gillespie received negative rxn_entry_mat")
    return false
  end
  if (T_sim<0)
    error("gillespie received negative T_sim")
    return false
  end
  return true
end

using Distributions

function gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim)
  if(!gillespie_input_ok(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim))
    error("gillespie received bad input.")
  else
    t_spent = 0
    rxn_times = Array(Float64, 0)
    rxn_types = Array(Int64, 0)
    current_x = init_x
    num_rxns_occ = 0
    num_molecule_types = size(sto_mat)[1]
    num_rxn_types = size(sto_mat)[2]

    while(t_spent < T_sim)
      #Get reaction propensities, prod_j c_i* (X_j choose k_ij)

      alpha = copy(rxn_rates)
      for rxn_index = [1:num_rxn_types]
        for mol_index = [1:num_molecule_types]
          alpha[rxn_index] = alpha[rxn_index]*binomial(current_x[mol_index], rxn_entry_mat[mol_index, rxn_index])
        end
      end
      alpha_sum = sum(alpha)

      #Perchance we're totally out of molecules, stop the simulation.
      #(It'll stop: see next block)
      #Otherwise, step forward in time
      if alpha_sum == 0
        t_spent = T_sim + 0.00000001 #Overshoot microscopically just in case I put a > instead of >=
      else
        tau = rand(Exponential(1/alpha_sum))
        t_spent = t_spent + tau
      end


      #If time's not up, update the molecule counts, the reactions
      #count, and the reaction times and types list
      if(t_spent < T_sim)
        current_rxn_type = rand(Categorical(alpha/alpha_sum))
        current_x = current_x + sto_mat[:,current_rxn_type]
        num_rxns_occ = num_rxns_occ + 1

        #! means that this modifies the array that it is given
        push!(rxn_times, t_spent)
        push!(rxn_types, current_rxn_type)
      end
    end
  end
  return current_x, num_rxns_occ, rxn_types,rxn_times, t_spent
end



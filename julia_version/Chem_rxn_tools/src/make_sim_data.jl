function make_sim_data(t_obs, true_cri, obs_mol_name, noise_distribution)
  T_sim = t_obs[end]
  #-----------------------------Simulate the complete data-------------------------------
  sim_results = Chem_rxn_tools.gillespie(true_cri, T_sim, false)
  sim_results.obs_mol_name = obs_mol_name
  sim_results.obs_mol_ind = Chem_rxn_tools.get_chem_indices(true_cri, obs_mol_name)

  #-----------------------------Simulate the observations-------------------------------
  sim_results.t_obs = t_obs
  x_index_matching_t_obs = Array(Int64, length(t_obs))
  sim_results.x_obs = Array{Int64,1}[]
  sim_results.d_obs = zeros(Float64,length(t_obs)) #Observed values
  for j in 1:length(t_obs)
    #Assume at first that t_obs[j] is after all rxns have happened.
    push!(sim_results.x_obs, sim_results.x_path[sim_results.num_rxns_occ])
    #Test other possible windows between reaction times. If
    #one of them works out, exit the loop knowing that t_obs[j] includes only first i-1 rxns
    for i in 2:sim_results.num_rxns_occ
      if(sim_results.rxn_times[i-1]<=t_obs[j] && t_obs[j] <= sim_results.rxn_times[i])
        x_index_matching_t_obs[j] = i - 1
        sim_results.x_obs[j] = sim_results.x_path[x_index_matching_t_obs[j]]
        break
      end
    end
    sim_results.d_obs[j] = sim_results.x_obs[j][sim_results.obs_mol_ind] + rand(noise_distribution, 1)[1]
  end
  return sim_results
end

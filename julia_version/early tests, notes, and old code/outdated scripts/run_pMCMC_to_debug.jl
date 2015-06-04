#This allows Julia to find Eric's modules.
workspace()
push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")

#-----------------------------Load in info from SBML shorthand and visualize results.------------------------------
  using Chem_rxn_tools
  SBML_file = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/wilkinson_rxns_SBML_shorthand.txt"
  wilk_cri = Chem_rxn_tools.SBML_read(SBML_file)
  sto_mat_graph,rxn_entry_mat_graph = Chem_rxn_tools.make_cri_graphic(wilk_cri)
  sto_mat_graph
  rxn_entry_mat_graph

#-----------------------------Set up the experiment-specific preferences-------------------------------
#Note: nothing in this code enforces that the same model is used to add noise and to calculate the logpdf of the noise.
#You can, purposefully or otherwise, run this code while simulating data from one model and doing inference on another.
#However, changing this line in this script as written May 17, 2015 will change both the simulated data and the inference method.
#The variable noise_distribution goes into both.
  using Distributions
  noise_distribution = Normal(0, 10)

#In Wilkinson's first experiment, the molecule of interest, SigD, is observed.
  obs_mol_name = "SigD"
  obs_mol_ind = Chem_rxn_tools.get_chem_indices(wilk_cri, obs_mol_name)

#In Wilkinson's first experiment, three rates are treated as unknown. Listed with ground truth values, they are:
#kSigDprod = 1, kflacherep = 0.02, kflacheunrep = 0.1.
  unk_names = ["SigDprod", "flacherep", "flacheunrep"]
  unk_rates = [1.0, 0.02, 0.1]
  unk_inds = Int64[]
  for word in unk_names
    push!(unk_inds, Chem_rxn_tools.get_rate_indices(wilk_cri, word))
  end
  #Sanity check:
  for i in 1:3
    println(string(unk_names[i], " with rate ", unk_rates[i], " should be: ", wilk_cri.rxn_labels[unk_inds[i]], " with rate ", wilk_cri.rxn_rates[unk_inds[i]] ,"."))
  end

  t_interval = 300.0
  num_intervals = 24.0
  T_sim = t_interval*num_intervals
  t_obs = t_interval*[1:num_intervals]

#-----------------------------Set up the sampler prefs and the prior-------------------------------

using pMCMC_julia
  MCS = pMCMC_julia.MCMC_state()
  MCS.bandwidth = 0.001
  MCS.do_kde = false
  num_samples_desired = 1000
  MCS.burnin_len = 1e3
  MCS.thin_len = 5
  param_sample = 10.^(-4*rand(3,num_samples_desired)) #log uniform prior on 1*10^-4 to 1*10^0 for the 3 unknowns
  state_sample = repmat(wilk_cri.init_amts,1,num_samples_desired)
  MCS.current_sample = pMCMC_julia.Sample_state_and_params_type(param_sample, state_sample)

  #Set up the model for emissions
  function emission_logden(x_current, d_obs)
    log_density = logpdf(noise_distribution, d_obs-x_current[obs_mol_ind])
    return log_density
  end
  MCS.emission_logden = emission_logden

  println("WARNING: fwd_sim WILL MODIFY wilk_cri WHEN pMCMC_julia.pMCMC IS CALLED.")
  println("It speeds things up not to copy the whole object.")
  println("Likely consequence: screwing up the init.amts field.")
  function fwd_sim(prev_sample_state::Array{Int64, 1}, prop_sample_params::Array{Float64, 1}, T_sim::Float64)
    for i in 1:length(unk_inds)
      wilk_cri.rxn_rates[unk_inds[i]] = prop_sample_params[i]
    end
    wilk_cri.init_amts = prev_sample_state
    return Chem_rxn_tools.gillespie(wilk_cri, T_sim, true)
  end
  MCS.fwd_sim = fwd_sim

#-----------------------------Save sampler prefs-------------------------------
using Dates
  metadata_to_save = string("This test was run at time ", now(), " with ",
                            "observations at intervals of ", t_interval,
                            " from time zero to ", T_sim,
                            " with a noise sd of ", std(noise_distribution), ". ",
                            " There were ",   MCS.current_sample.num_particles, " particles, ",
                            " with a log-uniform prior between 1 and 1e-4 and ",
                            " true rates of ", wilk_cri.rxn_rates, ". \n ",
                            "The observed molecule was ", obs_mol_name, " at index ", obs_mol_ind, ". \n"
                          )

#-----------------------------Simulate the hidden states-------------------------------
sim_results = Chem_rxn_tools.gillespie(wilk_cri, T_sim, false)

#-----------------------------Simulate the observations-------------------------------
x_index_matching_t_obs = Array(Int64, length(t_obs))
  x_obs = Int64[] #Noiseless verion of observed values.
  d_obs = zeros(Float64,length(t_obs)) #Observed values
  for j in 1:length(t_obs)
    #Assume at first that t_obs[j] is after all rxns have happened.
    push!(x_obs, sim_results.x_path[sim_results.num_rxns_occ][obs_mol_ind])
    #Test other possible windows between reaction times. If
    #one of them works out, exit the loop knowing that t_obs[j] includes only first i-1 rxns
    for i in 2:sim_results.num_rxns_occ
      if(sim_results.rxn_times[i-1]<=t_obs[j] && t_obs[j] <= sim_results.rxn_times[i])
        x_index_matching_t_obs[j] = i - 1
        x_obs[j] = sim_results.x_path[x_index_matching_t_obs[j]][obs_mol_ind]
        break
      end
    end
    d_obs[j] = x_obs[j] + rand(noise_distribution, 1)[1]
  end

#-----------------------------plot the simulated trajectory-------------------------------
  using Winston
  x_obs_vals = Array(Int64,length(t_obs))
  for i in 1:length(x_obs)
    x_obs_vals[i] = x_obs[i]
  end
  sim_trajec_plot = FramedPlot(
           title="Hidden (black) and observed (red) states from simulation",
           xlabel="Time",
           ylabel="State")
  add(sim_trajec_plot, Points(t_obs, x_obs_vals))
  add(sim_trajec_plot, Points(t_obs, d_obs, "color","red"))


#-----------------------------do the inference-------------------------------
using ProfileView
  are_we_profiling = false #faster not to profile
  tic()
  if(are_we_profiling)
    Profile.clear()
    @profile pMCMC_julia.pMCMC!(d_obs, t_obs, MCS)
    ProfileView.view()
  else
    pMCMC_julia.pMCMC!(d_obs, t_obs, MCS)
  end
  time_taken = toc()
  metadata_to_save = string(metadata_to_save, " It took ", time_taken, " seconds.")

today_filepath = string("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/experiments_after_tidy/debug_", now())
mkdir(today_filepath)

#plot the posterior
post_hists = pMCMC_julia.plot_from_MCMC_and_save(MCS, today_filepath, unk_rates, sim_results.x_path[end])

#-----------------------------save the results------------------------------
using HDF5, JLD
#To save important metadata
  save(string(today_filepath, "/metadata"), "metadata_to_save", metadata_to_save, "wilk_cri", wilk_cri)
  save(string(today_filepath, "/samples_and_data"), "posterior_sample", MCS.current_sample, "sim_results", sim_results, "x_obs", x_obs, "d_obs", d_obs, "t_obs", t_obs)



#This allows Julia to find Eric's modules.
workspace()
push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")
num_samples_desired = 1_000

#-----------------------------Load in info from SBML shorthand and visualize results.------------------------------
using Chem_rxn_tools
  demo_cri = Chem_rxn_info(
    "kryptonite",                        #species_labels
    50,                                  #init_amts
    1,                                   #num_species

    zeros(Int64, 1, 1),                  #sto_mat
    zeros(Int64, 1, 1),                  #rxn_entry_mat
    "Dummy_rxn",                         #rxn_labels
    1,                                   #num_rxns
    zeros(Float64,1),                    #rxn_rates
    "kryptonite->kryptonite",            #rxns_written_out

    [0],                                 #rxn_pos_in_SBML_file
    ["kryptonite_persistence"],          #SBML_par_names
    [1]                                  #rxn_labels
    )

  unk_inds = [1]
  unk_names = ["kryptonite_persistence"]
  obs_mol_name = "kryptonite"

  t_interval = 300.0
  num_intervals = 24.0
  T_sim = t_interval*num_intervals
  t_obs = t_interval*[1:num_intervals]
#-----------------------------Set up the sampler prefs and the prior-------------------------------

using pMCMC_julia
  MCS = pMCMC_julia.MCMC_state()
  MCS.bandwidth = 1
  MCS.do_kde = false
  MCS.burnin_len = 1e3
  MCS.thin_len = 5

  prior_sample_thetas = zeros(1,num_samples_desired)
  prior_sample = Sample_state_and_params_type(prior_sample_thetas, reshape(sample(1:500, num_samples_desired),1,num_samples_desired))
  MCS.current_sample = pMCMC_julia.Sample_state_and_params_type(param_sample, state_sample)

  obs_mol_ind = Chem_rxn_tools.get_chem_indices(demo_cri, obs_mol_name)

  using Distributions
  noise_distribution = Binomial(40, 0.5)
  #Set up the model for emissions
  function emission_logden(x_current, d_obs)
    log_density = logpdf(noise_distribution, d_obs-x_current[obs_mol_ind])
    return log_density
  end
  MCS.emission_logden = emission_logden

  function fwd_sim(prev_sample_state::Array{Int64, 1}, prop_sample_params::Array{Float64, 1}, T_sim::Float64)
    for i in 1:length(unk_inds)
      demo_cri.rxn_rates[unk_inds[i]] = prop_sample_params[i]
    end
    demo_cri.init_amts = prev_sample_state
    return Chem_rxn_tools.gillespie(demo_cri, T_sim, true)
  end
  MCS.fwd_sim = fwd_sim

#-----------------------------Save sampler prefs-------------------------------
using Dates
  metadata_to_save = string("This test was run at time ", now(), " with ",
                            "observations at intervals of ", t_interval,
                            " from time zero to ", T_sim,
                            " with a noise sd of ", std(noise_distribution), ". ",
                            " There were ",   MCS.current_sample.num_particles, " particles, ",
                        " with a clairvoyant spike-at-zero prior on reaction rate, ",
                        " a discrete unif(0,500) prior on the state, ",
                        " and true state of ", true_init_x, ".")

today_filepath = string("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/experiments_after_tidy/may27_verify", now())
  mkdir(today_filepath)

#-----------------------------Simulate the observations; plot; save plot-------------------------------
sim_results = Chem_rxn_tools.make_sim_data(t_obs, demo_cri, obs_mol_name, noise_distribution)
  mols_to_show = ["kryptonite"]
  Chem_rxn_tools.plot_save_sim_data(today_filepath, sim_results, demo_cri, mols_to_show)

#-----------------------------do the inference-------------------------------
using ProfileView
  are_we_profiling = false #faster not to profile
  tic()
  if(are_we_profiling)
    Profile.clear()
    @profile pMCMC_julia.pMCMC!(sim_results.d_obs, sim_results.t_obs, MCS)
    ProfileView.view()
  else
    pMCMC_julia.pMCMC!(sim_results.d_obs, sim_results.t_obs, MCS)
  end
  time_taken = toc()
  metadata_to_save = string(metadata_to_save, " It took ", time_taken, " seconds.")

#-----------------------------plot the posterior, print stuff-------------------------------
post_hists = pMCMC_julia.plot_save_marginals(MCS, today_filepath, unk_rates, sim_results.x_path[end])
  println(metadata_to_save)
  println(string("By round, the acceptance rate was ", MCS.num_acc/(MCS.burnin_len + MCS.thin_len*num_samples_desired)))

#-----------------------------save the results------------------------------
using HDF5, JLD
#To save important metadata
  save(string(today_filepath, "/metadata"), "metadata_to_save", metadata_to_save, "wilk_cri", wilk_cri)
  save(string(today_filepath, "/samples_and_data"), "num_acc", MCS.num_acc, "posterior_sample", MCS.current_sample, "sim_results", sim_results)




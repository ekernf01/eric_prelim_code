#This allows Julia to find Eric's modules.
workspace()
push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")

#The most oft changed preference
num_samples_desired = 10_000

#-----------------------------Load in info from SBML shorthand and visualize results.------------------------------
  using Chem_rxn_tools
  demo_cri = Chem_rxn_tools.make_demo_cri(2)
  sto_mat_graph,rxn_entry_mat_graph = Chem_rxn_tools.make_cri_graphic(demo_cri)
  sto_mat_graph
  rxn_entry_mat_graph

#-----------------------------Set up the experiment-specific preferences-------------------------------
#Note: nothing in this code enforces that the same model is used to add noise and to calculate the logpdf of the noise.
#You can, purposefully or otherwise, run this code while simulating data from one model and doing inference on another.
#However, changing this line in this script as written May 17, 2015 will change both the simulated data and the inference method.
#The variable noise_distribution goes into both.
using Distributions
  noise_distribution = Normal(0, 10)
  obs_mol_name = "mol1"
  num_unks = 2
  unk_names = ["prod_mol1", "decay_mol2"]
  unk_rates = [0.01, 0.0001]
  unk_inds = Int64[]
  for word in unk_names
    push!(unk_inds, Chem_rxn_tools.get_rate_indices(demo_cri, word))
  end
  #Sanity check:
  for i in 1:num_unks
    println(string(unk_names[i], " with rate ", unk_rates[i], " should be: ", demo_cri.rxn_labels[unk_inds[i]], " with rate ", demo_cri.rxn_rates[unk_inds[i]] ,"."))
  end

  t_interval = 300.0
  num_intervals = 24.0
  T_sim = t_interval*num_intervals
  t_obs = t_interval*[1:num_intervals]

#-----------------------------Set up the sampler prefs and the prior-------------------------------

using pMCMC_julia
  MCS = pMCMC_julia.MCMC_state()
  MCS.bandwidth_multiplier = 0.1
  MCS.do_kde = true
  MCS.burnin_len = 1e3
  MCS.thin_len = 5
  param_sample = zeros(num_unks,num_samples_desired)
  param_sample[1, :] = 10.^(-4*rand(1,num_samples_desired))
  param_sample[2, :] = 10.^(-2-4*rand(1,num_samples_desired))
  state_sample = repmat(demo_cri.init_amts,1,num_samples_desired)
  MCS.current_sample = pMCMC_julia.Sample_state_and_params_type(param_sample, state_sample)

  obs_mol_ind = Chem_rxn_tools.get_chem_indices(demo_cri, obs_mol_name)
  #Set up the model for emissions
  function emission_logden(x_current, d_obs)
    log_density = logpdf(noise_distribution, d_obs-x_current[obs_mol_ind])
    return log_density
  end
  MCS.emission_logden = emission_logden

  println("WARNING: fwd_sim WILL MODIFY demo_cri WHEN pMCMC_julia.pMCMC IS CALLED.")
  println("It speeds things up not to copy the whole object.")
  println("Likely consequence: screwing up the init.amts field.")
  function fwd_sim(prev_sample_state::Array{Int64, 1}, prop_sample_params::Array{Float64, 1}, T_sim::Float64)
    for i in 1:length(unk_inds)
      demo_cri.rxn_rates[unk_inds[i]] = prop_sample_params[i]
    end
    demo_cri.init_amts = prev_sample_state
    return Chem_rxn_tools.gillespie(demo_cri, T_sim, true)
  end
  MCS.fwd_sim = fwd_sim

#-----------------------------setup what and where to save-------------------------------
using Dates
  MCS.save_path = string("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/replication_exps/demo_system", now(), "/")
  mkdir(MCS.save_path)
  metadata_to_save = string("This test was run at time ", now(), " with ",
                            "observations at intervals of ", t_interval,
                            " from time zero to ", T_sim,
                            " with a noise sd of ", std(noise_distribution), ". ",
                            " There were ",   MCS.current_sample.num_particles, " particles, ",
                            " with log-uniform priors for the three unknowns and ",
                            " true rates of ", demo_cri.rxn_rates, ". \n ",
                            "The observed molecule was ", obs_mol_name, " at index ", Chem_rxn_tools.get_chem_indices(demo_cri, obs_mol_name), ". \n",
                            "The purpose of this test was to verify the pMCMC is correct. The bandwidth multiplier was ", MCS.bandwidth_multiplier, ". "
                            )



#-----------------------------Simulate the observations; plot; save plot-------------------------------
sim_results = Chem_rxn_tools.make_sim_data(t_obs, demo_cri, obs_mol_name, noise_distribution)
  mols_to_show = ["mol1", "mol2"]
  Chem_rxn_tools.plot_save_sim_data(MCS.save_path, sim_results, demo_cri, mols_to_show)

using HDF5, JLD
  @save string(MCS.save_path, "/everything_just_before_inference")
  pMCMC_julia.MCS_save("MCS_just_before_inference", MCS)

#-----------------------------do the inference-------------------------------
using ProfileView
  are_we_profiling = true #marginally faster not to profile
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

#-----------------------------save the results------------------------------
using HDF5, JLD
  pMCMC_julia.MCS_save("samples_and_metadata", MCS)

  println(metadata_to_save)
  println(string("By round, the acceptance rate was ", MCS.num_acc/(MCS.burnin_len + MCS.thin_len*num_samples_desired)))

#-----------------------------plot the posterior-------------------------------

@load string(MCS.save_path, "/samples_and_metadata")
  post_hists = pMCMC_julia.plot_save_marginals(MCS, MCS.save_path, unk_rates, sim_results.x_path[end])
  include("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/contour_bivariate_plot_maker.jl")
  post_contour = contour_plot_two_mols(MCS, MCS.save_path, unk_names[num_unks-1], unk_names[num_unks], demo_cri, unk_rates, unk_names)
  include("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/bivariate_plot_maker.jl")
  post_biv = plot_save_two_mols(MCS, MCS.save_path, unk_names[num_unks-1], unk_names[num_unks], demo_cri, unk_rates, unk_names)
  include("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/plot_biv_at_each_stage.jl")
  plot_biv_at_each_stage(MCS.save_path, unk_names[num_unks-1], unk_names[num_unks], unk_rates, unk_names)


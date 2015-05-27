#This allows Julia to find Eric's modules.
workspace()
push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")
num_samples_desired = 100_000

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
  MCS.bandwidth = 0.05
  MCS.do_kde = true
  MCS.burnin_len = 1e3
  MCS.thin_len = 5
#log uniform prior on 1*10^-4 to 1*10^0 for the 3 unknowns:
  #log(kSigDprod) ∼ Unif(log{0.01}, log{100}),
  #log(kflacherep) ∼ Unif(log{0.0002}, log{2}),
  #log(kflacheunrep) ∼ Unif(log{0.001}, log{10}).
  param_sample = zeros(3,num_samples_desired)
  param_sample[1, :] = 10.^(2 - 4*rand(1,num_samples_desired))
  param_sample[2, :] = 2*(10.^(-4*rand(1,num_samples_desired)))
  param_sample[3, :] = 10.^(1 - 4*rand(1,num_samples_desired))
  state_sample = repmat(wilk_cri.init_amts,1,num_samples_desired)
  MCS.current_sample = pMCMC_julia.Sample_state_and_params_type(param_sample, state_sample)

  obs_mol_ind = Chem_rxn_tools.get_chem_indices(wilk_cri, obs_mol_name)
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
                            "The observed molecule was ", obs_mol_name, " at index ", Chem_rxn_tools.get_chem_indices(wilk_cri, obs_mol_name), ". \n")

today_filepath = string("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/experiments_after_tidy/SigD_obs_", now())
  mkdir(today_filepath)

#-----------------------------Simulate the observations; plot; save plot-------------------------------
sim_results = Chem_rxn_tools.make_sim_data(t_obs, wilk_cri, obs_mol_name, noise_distribution)
  mols_to_show = ["SigD", "CodY", "Hag"]
  Chem_rxn_tools.plot_save_sim_data(today_filepath, sim_results, wilk_cri, mols_to_show)

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

#-----------------------------plot the posterior-------------------------------
post_hists = pMCMC_julia.plot_save_marginals(MCS, today_filepath, unk_rates, sim_results.x_path[end])
include("bivariate_plot_maker.jl")
post_biv = plot_save_two_mols(MCS, today_filepath, unk_names[2], unk_names[3], wilk_cri, unk_rates, unk_names)

println(metadata_to_save)
println(string("By round, the acceptance rate was ", MCS.num_acc/(MCS.burnin_len + MCS.thin_len*num_samples_desired)))

#-----------------------------save the results------------------------------
using HDF5, JLD
#To save important metadata
  save(string(today_filepath, "/metadata"), "metadata_to_save", metadata_to_save, "wilk_cri", wilk_cri)
  save(string(today_filepath, "/samples_and_data"), "num_acc", MCS.num_acc, "posterior_sample", MCS.current_sample, "sim_results", sim_results)




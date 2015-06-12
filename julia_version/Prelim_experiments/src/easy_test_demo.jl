workspace()
  push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")
  using Prelim_experiments

ep = Prelim_experiments.Exp_prefs(
    1_000_000,#num_samples_desired,
    2.5, #bw_max,
    0.001, #bw_min,
    0.1, #bandwidth_multiplier,
    0.001,     #true_noise_sd, #for bureaucratic reasons, cannot quite set this to zero
    10,      #assumed_noise_sd,
    "demo", #cri_source,
    "mol1", #obs_mol_name,
    ["mol1"],  #mols_to_show
    ["decay", "prod"], #unk_names,
    10, #t_interval,
    50, #num_intervals,
    "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/replication_exps/easy_tests" #save_folder
    )

  #use this line to see a printout of everything the sampler is doing.
  #ep.verbose = true

using ProfileView
  Profile.clear()
  @profile Prelim_experiments.start_off_test(ep, ep.num_intervals)
  ProfileView.view()

#Prelim_experiments.make_all_plots(joinpath(ep.save_folder, string(ep.time_of_test)))


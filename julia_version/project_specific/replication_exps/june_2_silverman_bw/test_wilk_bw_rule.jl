workspace()
  push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")
  using Prelim_experiments
  ep = Prelim_experiments.Exp_prefs(
    1_000_000,#num_samples_desired,
    2.5, #bw_max,
    0.001, #bw_min,
    0.1, #bandwidth_multiplier,
    10, #noise_sd,
    "wilk", #cri_source,
    "SigD", #obs_mol_name,
    ["SigDprod", "flacherep", "flacheunrep"], #unk_names,
    300.0, #t_interval,
    24, #num_intervals,
    "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/replication_exps/june_2_silverman_bw" #save_folder
    )

  using ProfileView
  Profile.clear()
  @profile Prelim_experiments.start_off_test(ep, 1)
  ProfileView.view()

#Prelim_experiments.make_all_plots(joinpath(ep.save_folder, string(ep.time_of_test)))


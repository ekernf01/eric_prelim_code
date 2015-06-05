workspace()
  push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")
  using Prelim_experiments
  ep = Prelim_experiments.Exp_prefs(
    10,#num_samples_desired,
    1, #bw_max,
    1, #bw_min,
    1, #bandwidth_multiplier,
    10, #noise_sd,
    "wilk", #cri_source,
    "SigD", #obs_mol_name,
    ["SigDprod", "flacherep", "flacheunrep"], #unk_names,
    300.0, #t_interval,
    6, #num_intervals,
    "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/Prelim_experiments/test_results" #save_folder
    )

  Prelim_experiments.run_test_generic(ep)



workspace()
  push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")
  using Prelim_experiments
  ep = Prelim_experiments.Exp_prefs(
    10,#num_samples_desired,
    1, #bw_max,
    1, #bw_min,
    0.1, #bandwidth_multiplier,
    10, #noise_sd,
    "wilk", #cri_source,
    "SigD", #obs_mol_name,
    ["SigDprod", "flacherep", "flacheunrep"], #unk_names,
    300.0, #t_interval,
    24, #num_intervals,
    "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/Prelim_experiments/src/tests" #save_folder
    )


using ProfileView
Profile.clear()
  @profile Prelim_experiments.start_off_test(ep, 1)

  ProfileView.view()

last_test_load_path = joinpath(ep.save_folder, string(ep.time_of_test))

using ProfileView
Profile.clear()
  @profile Prelim_experiments.pick_up_test(last_test_load_path, 1)
  ProfileView.view()


workspace()
  push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")
  using Prelim_experiments
  using Distributions
  using Winston
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

  noise_distribution = Normal(0, ep.noise_sd)
  T_sim = ep.t_interval*ep.num_intervals
  t_obs = ep.t_interval*[1:ep.num_intervals]

  SBML_file = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/wilkinson_rxns_SBML_shorthand.txt"
  wilk_cri = Chem_rxn_tools.SBML_read(SBML_file)

  #-----------------------------Simulate the observations; plot; save plot-------------------------------
  sim_results = Chem_rxn_tools.make_sim_data(t_obs, wilk_cri, ep.obs_mol_name, noise_distribution)

println(names(sim_results))
n = sim_results.num_rxns_occ/ep.num_intervals
proportions = [sum([sim_results.rxn_types[j]==i for j in 1:n]) for i in 1:18]
M = ([[sim_results.x_path[i][j] for i in 1:n] for j in 1:13])
[maximum(M[i]) for i in 1:13]
plot(M[1])


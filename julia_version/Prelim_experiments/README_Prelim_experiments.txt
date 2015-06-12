README_Prelim_experiments

Prelim_experiments is a module that helps manage the overall structure of the project. The main data type used here is Exp_prefs, with instances usually named ep. The easiest interfaces to use are the script easy_test_demo, which runs an experiment with a small demo system, and the easy_test_wilk script, which runs the first of Darren Wilkinson's experiments.

Documentation for data type Exp_prefs:
  num_samples_desired: how many particles to use in the pMCMC.
  bw_max: for choice of bandwidth in pMCMC. See pMCMC module for further documentation.
  bw_min: for choice of bandwidth in pMCMC. See pMCMC module for further documentation.
  bandwidth_multiplier: for choice of bandwidth in pMCMC. See pMCMC module for further documentation.

  true_noise_sd: noise sd to use when generating synthetic data.
  ass_noise_sd: noise sd to assume when running the inference.

  cri_source: either "wilk" for the Bacillus subtilis motility regulation network or "demo" for a single-molecule system.
  verbose: For the pMCMC sampler. See pMCMC module for further documentation.

  obs_mol_name: what molecule to observe.
  obs_mol_ind: its index, integer; internal member, set automatically by start_off_test
  mols_to_show: which molecules to plot when displaying the simulated data.

  unk_names: names of reaction rates to be treated as unknown.
  unk_rates: ground truth values for unknown reaction rates.
  unk_inds: internal member, set automatically by start_off_test
  num_unks: internal member, set automatically by start_off_test

  t_interval: length of intervals at which to gather data
  num_intervals: number of data points gathered
  save_folder: path to save results and plots to

  time_of_test: internal member, set automatically by start_off_test
  time_taken: internal member, set automatically by start_off_test



start_off_test(ep, num_stages_to_try)
pick_up_test(save_path, num_stages_to_try)
These tests can take a while. You may want to start off a test, but tell it to stop after a few stages. This pair of functions will let you start it off, then pick it up again. This function will start at the last stage you did, and it will do num_stages_to_try more of them. The pMCMC_julia module has no native ability to do only the first few stages. You have to trick it: feed it only a segment of the data at first. start_off_test() does this automatically. The rest of the information about these functions appears in the docs for the Exp_prefs data type.

make_all_plots(save_path)--given a folder where start_off_test() (and maybe pick_up_test()) has put results, make_all_plots will sift through the saved samples from intermediate and final stages, re-doing heatmaps of the posterior densities.
plot_biv_at_each_stage(today_filepath)--a subroutine of make_all_plots().

plot_save_two_mols--deprecated. Use at our own risk.
contour_plot_two_mols--deprecated. Use at our own risk.
run_test_generic--deprecated. Use at our own risk.
plot_save_trajec_CIs--deprecated. Use at our own risk.
plot_save_marginals--deprecated. Use at our own risk.

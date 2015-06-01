using HDF5, JLD
pwd()
today_filepath = "Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/replication_exps/SigD_obs_2015-05-30T15:54:57"
cd(today_filepath)
readdir()

load("samples_and_metadata")
  post_hists = pMCMC_julia.plot_save_marginals(MCS, today_filepath, unk_rates, sim_results.x_path[end])
  include("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/contour_bivariate_plot_maker.jl")
  post_contour = contour_plot_two_mols(MCS, today_filepath, unk_names[2], unk_names[3], wilk_cri, unk_rates, unk_names)
  include("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/bivariate_plot_maker.jl")
  post_biv = plot_save_two_mols(MCS, today_filepath, unk_names[2], unk_names[3], wilk_cri, unk_rates, unk_names)


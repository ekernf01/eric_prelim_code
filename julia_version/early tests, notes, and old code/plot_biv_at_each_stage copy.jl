using Chem_rxn_tools
using pMCMC_julia
using HDF5, JLD
using Winston
using KernelDensity

function plot_biv_at_each_stage(today_filepath::String)
  cd()
  cd(today_filepath)

  unk_names = load(string("everything_just_before_inference"), "unk_names")
  num_unks = load(string("everything_just_before_inference"), "num_unks")
  unk_rates = load(string("everything_just_before_inference"), "unk_rates")
  wilk_cri = load(string("everything_just_before_inference"), "wilk_cri")
  ts_len = length(load(string("everything_just_before_inference"), "t_obs"))
  rate_name_x = unk_names[end-1]
  rate_name_y = unk_names[end]
  folder_w_plots = joinpath(today_filepath, "stagewise_plots/")
  if !isdir(folder_w_plots)
    mkdir(folder_w_plots)
  end

  MCS = pMCMC_julia.MCMC_state()
  for i in 0:(ts_len-1)
    MCS.current_sample = load(string("stage", i, "sample"),"current_sample")

    rate_ind_x = findin(unk_names,[rate_name_x])[1]
    rate_ind_y = findin(unk_names,[rate_name_y])[1]
    x = log(MCS.current_sample.params[rate_ind_x, :][:])
    y = log(MCS.current_sample.params[rate_ind_y, :][:])
    posterior_plot = imagesc(kde((x, y)))

    title("Bivariate MCMC output, truth in red. Log(rate) displayed.")
    xlabel(rate_name_x)
    ylabel(rate_name_y)

    add(posterior_plot, Points(log(unk_rates[rate_ind_x]), log(unk_rates[rate_ind_y]), "color", "red"))
    savefig(string(today_filepath, "/stagewise_plots/dist",i,"_contour_", rate_name_x, "_", rate_name_y, ".png"))
  end
end


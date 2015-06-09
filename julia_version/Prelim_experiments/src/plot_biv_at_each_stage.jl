function plot_biv_at_each_stage(today_filepath::String)
  cd()
  cd(today_filepath)
  ep = load(string("everything_just_before_inference"), "ep")
  wilk_cri = load(string("everything_just_before_inference"), "wilk_cri")
  sim_results = load(string("everything_just_before_inference"), "sim_results")
  MCS = pMCMC_julia.MCS_load("after_inf_sample")

  rate_name_x = ep.unk_names[end-1]
  rate_name_y = ep.unk_names[end]
  folder_w_plots = joinpath(today_filepath, "stagewise_plots/")
  if !isdir(folder_w_plots)
    mkdir(folder_w_plots)
  end

  for i in 0:ep.num_intervals
#     println("plot_biv_at_each_stage making plot: ", i)
    if i==ep.num_intervals
      MCS = pMCMC_julia.MCS_load(string("after_inf_sample"))
    else
      MCS = pMCMC_julia.MCS_load(string("stage", i, "sample"))
    end
    rate_ind_x = findin(ep.unk_names,[rate_name_x])[1]
    rate_ind_y = findin(ep.unk_names,[rate_name_y])[1]
    x = log(MCS.current_sample.params[rate_ind_x, :][:])
    y = log(MCS.current_sample.params[rate_ind_y, :][:])
    gc()
    my_kde = kde((x, y))
    posterior_plot = imagesc(my_kde)
    gc()
    title("Bivariate MCMC output, truth in red. Log(rate) displayed.")
    xlabel(rate_name_x)
    ylabel(rate_name_y)

    add(posterior_plot, Points(log(ep.unk_rates[rate_ind_x]), log(ep.unk_rates[rate_ind_y]), "color", "red"))
    savefig(string(today_filepath, "/stagewise_plots/dist",i,"_contour_", rate_name_x, "_", rate_name_y, ".png"))
  end
end


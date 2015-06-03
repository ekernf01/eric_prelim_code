using Chem_rxn_tools
using pMCMC_julia
using HDF5, JLD
using Winston
using KernelDensity
#today_filepath = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/replication_exps/may29 large bandwidth tests/wasted_run"

function plot_biv_at_each_stage(today_filepath::String, rate_name_x, rate_name_y, ground_truth_params=[], unk_names=[])
  cd()
  cd(today_filepath)
  mkdir(joinpath(today_filepath, "stagewise_plots/"))

  MCS = pMCMC_julia.MCMC_state()
  for i in 1:23
    my_dict = load(string("stage", i, "sample"))
    MCS.current_sample = my_dict["current_sample"]

    rate_ind_x = findin(unk_names,[rate_name_x])[1]
    rate_ind_y = findin(unk_names,[rate_name_y])[1]
    x = log(MCS.current_sample.params[rate_ind_x, :][:])
    y = log(MCS.current_sample.params[rate_ind_y, :][:])
    posterior_plot = imagesc(kde((x, y), boundary=((-10,5), (-10,5))))

    title("Bivariate MCMC output, truth in red. Log(rate) displayed.")
    xlabel(rate_name_x)
    ylabel(rate_name_y)

    add(posterior_plot, Points(log(unk_rates[rate_ind_x]), log(unk_rates[rate_ind_y]), "color", "red"))
    savefig(string(today_filepath, "/stagewise_plots/dist",i,"_contour_", rate_name_x, "_", rate_name_y, ".png"))
  end
end

using Chem_rxn_tools
using pMCMC_julia
using Winston
using KernelDensity



function contour_plot_two_mols(MCS::pMCMC_julia.MCMC_state, save_folder::String, rate_name_x, rate_name_y, cri::Chem_rxn_tools.Chem_rxn_info,ground_truth_params=[], unk_names=[])
  rate_ind_x = findin(unk_names,[rate_name_x])[1]
  rate_ind_y = findin(unk_names,[rate_name_y])[1]
  x = log(MCS.current_sample.params[rate_ind_x, :][:])
  y = log(MCS.current_sample.params[rate_ind_y, :][:])
  posterior_plot = imagesc(kde((x, y)))
  title("Bivariate MCMC output, truth in red. Log(rate) displayed.")
  xlabel(rate_name_x)
  ylabel(rate_name_y)

  if !isempty(ground_truth_params)
    add(posterior_plot, Points(log(ground_truth_params[rate_ind_x]), log(ground_truth_params[rate_ind_y]), "color", "red"))
  end
  savefig(string(save_folder, "/post_contour_", rate_name_x, "_", rate_name_y, ".png"))
  return posterior_plot
end



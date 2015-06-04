using Chem_rxn_tools
using pMCMC_julia
using Winston
println("Warning: plot_save_two_mols in bivariate_plot_maker has been deprecated.")
function plot_save_two_mols(MCS::pMCMC_julia.MCMC_state, save_folder::String, rate_name_x, rate_name_y, cri::Chem_rxn_tools.Chem_rxn_info,ground_truth_params=[], unk_names=[])
  num_particles = MCS.current_sample.num_particles
  skip_len = ceil(Int64, num_particles/1000)
  plot_from_indices = 1:skip_len:num_particles
  rate_ind_x = findin(unk_names,[rate_name_x])[1]
  rate_ind_y = findin(unk_names,[rate_name_y])[1]
  x = log10(MCS.current_sample.params[rate_ind_x, plot_from_indices][:])
  y = log10(MCS.current_sample.params[rate_ind_y, plot_from_indices][:])
  posterior_plot = plot(x, y, "k^")
  title("Bivariate MCMC output, truth in red. Log_{10}(rate) displayed.")
  xlabel(rate_name_x)
  ylabel(rate_name_y)

  if !isempty(ground_truth_params)
    add(posterior_plot, Points(log10(ground_truth_params[rate_ind_x]), log10(ground_truth_params[rate_ind_y]), "color", "red"))
  end
  savefig(string(save_folder, "/post_biv_", rate_name_x, "_", rate_name_y, ".png"))
  return posterior_plot
end


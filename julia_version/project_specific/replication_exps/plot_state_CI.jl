#Script to plot true hidden state along with predictive intervals
using Chem_rxn_tools
using pMCMC_julia
using HDF5, JLD
using Winston
function plot_save_trajec_CIs(save_path, sim_results::Chem_sim_result, cri::Chem_rxn_info, mols_to_show)

  save_path = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/replication_exps/debugging/SigD_obs_2015-06-01T19:58:04"
  @load string(save_path,"/", "everything_just_before_inference")

  CI_trajec_plot = plot(Int64[0],Int64[0],
           title="Trajectories of several molecules with predictive intervals",
           xlabel="Time",
           ylabel="Count")
  x_to_display = zeros(Int64,length(sim_results.t_obs), length(mols_to_show))
  all_curves = Points[]
  color_list = ["black", "red", "blue", "green"]

  #Outer loop: make the plot, adding one molecule at a time
  for j in 1:length(mols_to_show)
    mol_name = mols_to_show[j]
    mol_ind = Chem_rxn_tools.get_chem_indices(cri, mol_name)
    #inner loop: move along the observations
    for i in 1:length(sim_results.t_obs)
      x_to_display[i,j] = sim_results.x_obs[i][mol_ind]
      @load string(save_path,"/", "stage", i, "sample")
      mol_count_post_sample = current_sample.state[mol_ind,:][:]
      ci_upper_bound = quantile(mol_count_post_sample, 0.95)
      ci_lower_bound = quantile(mol_count_post_sample, 0.05)
      ci_center[i,j] = 0.5*(ci_upper_bound + ci_lower_bound)
      ci_radius[i,j] = 0.5*(ci_upper_bound - ci_lower_bound)
    end
    push!(all_curves, Points(sim_results.t_obs, x_to_display[:, j], "color", color_list[j]))
    setattr(all_curves[j], "label", mol_name)
    add(CI_trajec_plot, all_curves[j])
    add(CI_trajec_plot, SymmetricErrorBarsY(t_obs, ci_center[:, j][:], ci_radius[:, j][:]))
  end
  MyLegend = Legend(0.1, 0.9,  all_curves)
  add(sim_trajec_plot, all_curves..., MyLegend)
  savefig(string(save_path, "/true_data_with_CIs.png"))
end


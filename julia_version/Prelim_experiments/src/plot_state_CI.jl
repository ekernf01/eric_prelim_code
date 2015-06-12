#Script to plot true hidden state along with predictive intervals
using Chem_rxn_tools
using pMCMC_julia
using HDF5, JLD
using Winston
function plot_save_trajec_CIs(save_path)

  file_with_metadata = joinpath(save_path, "everything_just_before_inference")
  loaded_dict = load(file_with_metadata)
  ep = loaded_dict["ep"]
  sim_results = loaded_dict["sim_results"]
  cri = loaded_dict["wilk_cri"]
  mols_to_show = ep.mols_to_show

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
    for i in 0:(length(sim_results.t_obs)-1)
      x_to_display[i,j] = sim_results.x_obs[i][mol_ind]
      loaded_dict_i = load(string(save_path,"/", "stage", i, "sample"))
      current_sample = loaded_dict_i["current_sample"]
      mol_count_post_sample = current_sample.state[mol_ind,:][:]
      current_sample = []
      gc()
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


function plot_save_sim_data(today_filepath, sim_results::Chem_sim_result, cri::Chem_rxn_info, mols_to_show)
  #Observed data and true value at that coordinate
  x_obs_coord = Array(Int64,length(sim_results.t_obs))
  for i in 1:length(sim_results.t_obs)
    x_obs_coord[i] = sim_results.x_obs[i][sim_results.obs_mol_ind]
  end
  sim_trajec_plot = plot(Int64[0],Int64[0],
           title="Observed data (red) and true value (black) at that coordinate",
           xlabel="Time",
           ylabel="State")
  add(sim_trajec_plot, Points(sim_results.t_obs, x_obs_coord))
  add(sim_trajec_plot, Points(sim_results.t_obs, sim_results.d_obs, "color","red"))
  savefig(string(today_filepath, "/simulated_data_obs.png"))

  #Trajectories of several molecules
  sim_trajec_plot = plot(Int64[0],Int64[0],
           title="Trajectories of several molecules",
           xlabel="Time",
           ylabel="Count")
  x_to_display = zeros(Int64,length(sim_results.t_obs), length(mols_to_show))
  all_curves = Points[]
  color_list = ["black", "red", "blue", "green"]
  for j in 1:length(mols_to_show)
    mol_name = mols_to_show[j]
    for i in 1:length(sim_results.t_obs)
      x_to_display[i,j] = sim_results.x_obs[i][Chem_rxn_tools.get_chem_indices(cri, mol_name)]
    end
    push!(all_curves, Points(sim_results.t_obs, x_to_display[:, j], "color", color_list[j]))
    setattr(all_curves[j], "label", mol_name)
    sim_trajec_plot = add(sim_trajec_plot, all_curves[j])
  end
  MyLegend = Legend(0.1, 0.9,  all_curves)
  add(sim_trajec_plot, all_curves..., MyLegend)
  savefig(string(today_filepath, "/simulated_data_multiple_mols.png"))
  return sim_trajec_plot
end


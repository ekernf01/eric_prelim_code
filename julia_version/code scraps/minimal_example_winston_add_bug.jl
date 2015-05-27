#Trajectories of several molecules
sim_trajec_plot = plot(Int64[0],Int64[0],
           title="Trajectories of several molecules",
           xlabel="Time",
           ylabel="Count")
all_curves = Points[]
mols_to_show = ["a", "b"]
t_obs = [1:24]
x_obs = reshape(1:48, 24, 2)
for j in 1:length(mols_to_show)
  mol_name = mols_to_show[j]
  push!(all_curves, Points(t_obs, x_obs[:, j]))
  setattr(all_curves[j], "label", mol_name);
  #sim_trajec_plot = add(sim_trajec_plot, all_curves[j])
end
MyLegend = Legend(0.5, 0.9,  {all_curves[1], all_curves[2]})
add(sim_trajec_plot, all_curves[1],all_curves[2], MyLegend);


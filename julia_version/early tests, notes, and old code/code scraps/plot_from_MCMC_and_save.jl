
function plot_from_MCMC_and_save(MCS::pMCMC_julia.MCMC_state, save_folder::String, ground_truth_params=[], ground_truth_state=[])
  num_particles = MCS.current_sample.num_particles
  par_dim = MCS.current_sample.par_dim
  state_dim = MCS.current_sample.state_dim

  skip_len = ceil(Int64, num_particles/1000)
  plot_from_indices = 1:skip_len:num_particles

  for i in 1:par_dim
    to_plot = MCS.current_sample.params[i,plot_from_indices]
    posterior_plot = plot(Int64[], Int64[], title=string("Marginal", i, "of MCMC output (pars)"))
    add(posterior_plot, Histogram(hist(vec(to_plot))...))
    if !isempty(ground_truth_params)
      #Add a line of the right height
      line_height = maximum(hist(vec(to_plot))[2])
      add(posterior_plot, Curve([ground_truth_params[i],ground_truth_params[i]], [line_height,0], "color", "red"))
    end
    savefig(string(save_folder, "/post_par_", i, ".png"))
  end
  for i in 1:state_dim
    to_plot = MCS.current_sample.state[i,plot_from_indices]
    posterior_plot = plot(Int64[], Int64[], title=string("Marginal", i, "of MCMC output (state)"))
    add(posterior_plot, Histogram(hist(vec(to_plot))...))
    if !isempty(ground_truth_state)
      #Add a line of the right height
      line_height = maximum(hist(vec(to_plot))[2])
      add(posterior_plot, Curve([ground_truth_state[i],ground_truth_state[i]], [line_height,0], "color", "red"))
    end
    savefig(string(save_folder, "/post_state_", i, ".png"))
  end
end


#This subroutine converts a large sample from P(params, state at time i-1|data to time i-1) into
#a large sample from P(params, state at time i|data to time i)
function pMCMC_single_stage!(samples_post_prev_stage,
                             num_samples_desired,
                             burnin_len,
                             thin_len,
                             d_obs_this_stage,#algo only needs the next data point, not all at once
                             noise_sd,
                             obs_molecule_index,
                             sto_mat,
                             rxn_entry_mat,
                             T_sim_this_stage, #needs to know how far to simulate
                             bandwidth)

  par_dim = size(samples_post_prev_stage.params)[1]
  state_dim = size(samples_post_prev_stage.state)[1]
  #new empty array
  samples_post_current_stage = Sample_state_and_params_type(num_samples_desired,par_dim, state_dim)

  #Loop setup
  #this is an independence proposal, so I could draw all the proposals and unifs
  #beforehand, but I think it's actually faster in Julia to do it in a loop.
  #initial sample to default to when rejecting
  num_particles = size(samples_post_prev_stage.params)[2]
  particle_index = sample([1:num_particles])
  prev_sample_params = samples_post_prev_stage.params[:,particle_index]
  particle_index
  prev_sample_state = samples_post_prev_stage.state[:,particle_index]
  chain_len = burnin_len + thin_len*num_samples_desired
  record_index = 0
  num_acc = 0

  prop_sample_params = Array(Float64, par_dim)
  prop_sample_state = Array(Int64, state_dim)
  range_num_part = [1:num_particles]
  kde_kernel = Normal(0,bandwidth)
  for chainstep in 1:chain_len
    #param proposal from previous stage
    prop_particle_index = sample(range_num_part)
    prop_sample_params = samples_post_prev_stage.params[:,particle_index]*exp(rand(kde_kernel))#do kde in log space
    prop_sample_state = samples_post_prev_stage.state[:,particle_index]

    #add perturbation as if sampling from a KDE
    prop_sample_params = prop_sample_params

    #state proposal from fwd simulation
    inside_sampler =true
    prop_sample_state = gillespie(prop_sample_state, sto_mat, rxn_entry_mat, prop_sample_params, T_sim_this_stage, inside_sampler)

    #accept if A > 1 or A > unif, i.e. log>0 or log>log(unif)
    log_acc_rat =
      log_measurement_density(prop_sample_state, noise_sd, obs_molecule_index, d_obs_this_stage) -
      log_measurement_density(prev_sample_state, noise_sd, obs_molecule_index, d_obs_this_stage)
    if (log_acc_rat > 0) || (log_acc_rat > log(rand(1)[1]))
      acc_sample_state = prop_sample_state
      acc_sample_params = prop_sample_params
      num_acc = num_acc + 1
    else
      acc_sample_state = prev_sample_state
      acc_sample_params = prev_sample_params
    end

    #record after burnin unless thinned out
    if ((chainstep>burnin_len) && (chainstep%thin_len==1))
      record_index = record_index + 1
      samples_post_current_stage.params[:,record_index] = acc_sample_params
      samples_post_current_stage.state[:,record_index] = acc_sample_state
    end
    prev_sample_state = acc_sample_state
    prev_sample_params = acc_sample_params
  end

  return samples_post_current_stage, num_acc
end

#This function loops over the dataset. At iteration i, it calls an MCMC-based
#subroutine to convert a large sample from P(params|data to time i-1) into
#a large sample from P(params|data to time i)
function pMCMC(d_obs, t_obs, prior_sample,
                             num_samples_desired,
                             burnin_len,
                             thin_len,
                             noise_sd,
                             obs_molecule_index,
                             sto_mat,
                             rxn_entry_mat,
                             bandwidth)

  samples_post_current_stage = Sample_state_and_params_type(prior_sample.params, prior_sample.state)
  I = length(d_obs)
  num_acc = zeros(Int64, I)
  for stage = 1:I #by stage, I mean how much data has been conditioned upon. At stage 2, we've conditioned on 2 data points.

    #fold in more data
    if stage>1
      T_sim_this_stage = t_obs[stage] - t_obs[stage]
    else
      T_sim_this_stage = t_obs[stage]
    end
    d_obs_this_stage = d_obs[stage]

    #re-run the sampler
    samples_post_current_stage,num_acc[stage] = pMCMC_single_stage!(samples_post_current_stage,
                             num_samples_desired,
                             burnin_len,
                             thin_len,
                             d_obs_this_stage,#algo only needs the next data point, not all at once
                             noise_sd,
                             obs_molecule_index,
                             sto_mat,
                             rxn_entry_mat,
                             T_sim_this_stage, #needs to know how far to simulate
                             bandwidth)
  end
  return samples_post_current_stage,num_acc
end

function plot_from_MCMC(posterior_sample, i)
  margin = abs(minimum(posterior_sample.params[2,:]))*0.1
  posterior_plot = FramedPlot(
            title="Bivariate output of MCMC sampler, rates 1 and 2",
            xlabel="Rate 1",
            ylabel="Rate 2",
            xrange=(minimum(posterior_sample.params[1,:])-margin, maximum(posterior_sample.params[1,:])+margin),
            yrange=(minimum(posterior_sample.params[2,:])-margin, maximum(posterior_sample.params[2,:])+margin)
            )
  num_particles = size(posterior_sample.params)[2]
  if num_particles > 1000
    skip_len = floor(Int64, num_particles/1000)
    add(posterior_plot, Points(posterior_sample.params[1,[1:skip_len:num_particles]],posterior_sample.params[2,:]))
  else
    add(posterior_plot, Points(posterior_sample.params[1,:],posterior_sample.params[2,:]))
  end
  return posterior_plot
end

[1:2:10]




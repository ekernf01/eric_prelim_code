#This subroutine converts a large sample from P(params, state at time i-1|data to time i-1) into
#a large sample from P(params, state at time i|data to time i)
function pMCMC_single_stage!(samples_post_current_stage,

                             num_samples_desired,
                             burnin_len,
                             thin_len,
                             d_obs_this_stage,#algo only needs the next data point, not all at once
                             obs_molecule_index,
                             sto_mat,
                             rxn_entry_mat,
                             T_sim_this_stage #needs to know how far to simulate
                             )

  samples_post_prev_stage = samples_post_current_stage
  samples_post_current_stage = 0 #don't want two references to the same array

  #Loop setup
  #this is an independence proposal, so I could draw all the proposals and unifs
  #beforehand, but I think it's actually faster in Julia to do it in a loop.
  prev_sample = sample(samples_post_prev_stage) #initial sample to default to when rejecting
  chain_len = burnin_len + 1 + (thin_len-1)*num_samples_desired

  for chainstep in 1:chain_len
    #param proposal from previous stage
    prop_sample = sample(sample_post_prev_stage)

    #state proposal from fwd simulation
    prop_sample.state = gillespie(prop_sample.state, sto_mat, rxn_entry_mat, rxn_rates, T_sim_this_stage, inside_sampler=true)

    #accept if A > 1 or A > unif, i.e. log>0 or log>log(unif)
    log_acc_rat =
      log_measurement_density(prop_sample.state, noise_sd, molecule_index=obs_molecule_index, d_obs=d_obs_this_stage) -
      log_measurement_density(prev_sample.state, noise_sd, molecule_index=obs_molecule_index, d_obs=d_obs_this_stage)
    if (log_acc_rat > 0) || (log_acc_rat > log(rand(1)))
      acc_sample = prop_sample
    else
      acc_sample = prev_sample
    end

    #record after burnin unless thinned out
    if ((chainstep>burnin_len) && (chainstep%thin_len==1))
      samples_post_current_stage[(chainstep+4)/5] = acc_sample
    end
    prev_sample = acc_sample
  end

  return samples_post_current_stage
end


#This function loops over the dataset. At iteration i, it calls an MCMC-based
#subroutine to convert a large sample from P(params|data to time i-1) into
#a large sample from P(params|data to time i)
function pMCMC(d_obs, t_obs, prior_sample,
                             num_samples_desired,
                             burnin_len,
                             thin_len,
                             obs_molecule_index,
                             sto_mat,
                             rxn_entry_mat)

  sample_post_current_stage = prior_sample
  I = length(d_obs)
  for stage = 1:I #by stage, I mean how much data has been conditioned upon. At stage 2, we've conditioned on 2 data points.

    #fold in more data
    if stage>1
      T_sim_this_stage = t_obs[stage] - t_obs[stage]
    else
      T_sim_this_stage = t_obs[stage]
    end
    d_obs_this_stage = d_obs[stage]

    #re-run the sampler
    sample_post_current_stage = pMCMC_single_stage!(samples_post_current_stage,
                             num_samples_desired,
                             burnin_len,
                             thin_len,
                             d_obs_this_stage,#algo only needs the next data point, not all at once
                             obs_molecule_index,
                             sto_mat,
                             rxn_entry_mat,
                             T_sim_this_stage #needs to know how far to simulate
                             )
  end
  return sample_post_current_stage
end

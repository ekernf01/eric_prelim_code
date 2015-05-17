# This module implements the particle MCMC method from Darren Wilkinson's 2009 paper
# "Parameter inference for stochastic kinetic models of bacterial gene regulation: a Bayesian approach to systems biology"
# The method is for Bayesian inference on hidden markov models of continuous time-series data.
# This code currently supports only state spaces consisting of ordered lists of integers and parameter spaces contained in R^n.
# The main call is like this:
#
# using pMCMC_julia
# pMCMC(d_obs, t_obs, current_MCMC_state::MCMC_state, emission_logden::Function, fwd_sim::Function)
#
# It requires as input:
# Data: a vector of observations d_obs and times t_obs when they were made
# Preferences for the sampler in the form of current_MCMC_state
# Samples from your prior distribution on the initial state and on the params in the form of current_MCMC_state's current_sample field
# A fast routine fwd_sim() to simulate your model of choice forwards for an arbitrary amount of (simlated) time.
#   It should take:
#     parameters, (latent) initial state, desired duration
#   It should return:
#     Hidden state after simulation of desired duration
# A fast routine emission_logden() to compute the log-density of the measurement error
#   It should take:
#     Hidden state x_ob that may have happened at a given time, observation/emission from the same time point d_ob
#   It should return:
#     log-density of the measurement error assuming x_ob gave rise to d_ob, or in other language, log of the emission density e(d_ob|x_ob)

#This module also contains functions to plot and save results.


module pMCMC_julia
#---------------------------------Sample_state_and_params_type--------------------------------------------
#Holds a big array of samples on parameters and hidden states
#Methods to construct, copy, and check desired invariants
immutable Sample_state_and_params_type
  params::Array{Float64,2}
  state::Array{Int64,2}
  par_dim::Int64
  state_dim::Int64
  num_particles::Int64

  function Sample_state_and_params_type(params::Array{Float64, 2}, state::Array{Int64, 2})
    new_object = Sample_state_and_params_type(params,
                                   state,
                                   size(params)[1], #par_dim
                                   size(state)[1],  #state_dim
                                   max(size(params)[2], size(state)[2]) #num_particles
                                   )
    Sample_state_and_params_type_data_check(new_object)
    return new_object
  end

  function Sample_state_and_params_type(num_particles::Int64, par_dim::Int64, state_dim::Int64)
    new_object = Sample_state_and_params_type(zeros(Float64, par_dim, num_particles),
                                   zeros(Int64, state_dim, num_particles),
                                   par_dim,
                                   state_dim,
                                   num_particles)
    Sample_state_and_params_type_data_check(new_object)
    return new_object
  end
end
    function Sample_state_and_params_type_data_check(Sample_state_and_params_type)
      mismatch1 = (par_dim != size(samples_post_prev_stage.params)[1])
      mismatch2 = (state_dim != size(samples_post_prev_stage.state)[1])
      mismatch3 = (num_particles != size(samples_post_prev_stage.params)[2])
      mismatch4 = (num_particles != size(samples_post_prev_stage.state)[2])
      if(mismatch1 | mismatch2 | mismatch3 | mismatch4)
        error("dimension data don't match in a Sample_state_and_params_type object.")
      end
      return true
    end

    function copy_sample(current_sample::Sample_state_and_params_type)
      return Sample_state_and_params_type(current_sample.params,
                                          current_sample.state,
                                          current_sample.par_dim,
                                          current_sample.state_dim,
                                          current_sample.num_particles)
    end

#---------------------------------MCMC_state--------------------------------------------
#Holds samples and preferences for the pMCMC algo
#Method to construct, check invariants, and run a single stage of the sampler
type MCMC_state
  num_samples_desired::Int64
  burnin_len::Int64
  thin_len::Int64
  bandwidth::Float64
  do_kde::Bool
  num_acc::Array{Int64, 1}
  stage::Int64
  current_sample::Sample_state_and_params_type
end
    MCMC_state() = MCMC_state(
      1000,    #num_samples_desired
      1000,    #burnin_len
      5,       #thin_len
      0.001,   #bandwidth
      false,   #do_kde
      Int64[], #num_acc
      0,       #stage
      Sample_state_and_params_type(0, 0, 0)
      )
    end

    function MCMC_state_data_check(current_sampler_state::MCMC_state)
      if length(current_sampler_state.num_acc!=stage)
        error("MCMC_state has mismatched fields: stage and num_acc")
      end
      Sample_state_and_params_type_data_check(current_sampler_state.current_sample)
      return true
    end
#This function loops over the dataset. At iteration i, it calls an MCMC-based
#subroutine to convert a large sample from P(params|data to time i-1) into
#a large sample from P(params|data to time i)
function pMCMC(d_obs, t_obs,
               current_MCMC_state::MCMC_state,
               emission_logden::Function,
               fwd_sim::Function)

  samples_post_current_stage = copy_sample(current_MCMC_state.current_sample)
  I = length(d_obs)

  for stage = 1:I #by stage, I mean how much data has been conditioned upon. At stage 2, we've conditioned on 2 data points.
    println(string("In pMCMC at stage ", stage, "."))
    Sample_state_and_params_type_data_check(current_MCMC_state.current_sample::Sample_state_and_params_type)
    #fold in more data
    if stage>1
      T_sim_this_stage = t_obs[stage] - t_obs[stage]
    else
      T_sim_this_stage = t_obs[stage]
    end
    d_obs_this_stage = d_obs[stage]

    #re-run the sampler
    pMCMC_single_stage!(current_MCMC_state::MCMC_state,
                        emission_logden::Function,
                        fwd_sim::Function)
  end
  return current_MCMC_state
end

#This subroutine converts a large sample from P(params, state at time i-1|data to time i-1) into
#a large sample from P(params, state at time i|data to time i)
function pMCMC_single_stage!(current_MCMC_state::MCMC_state,
                              emission_logden::Function,
                              fwd_sim::Function)

  #new empty array
  samples_post_current_stage = Sample_state_and_params_type(num_samples_desired,par_dim, state_dim)

  #Loop setup
  #this is an independence proposal, so I could draw all the proposals and unifs
  #beforehand, but I think it's actually faster in Julia to do it in a loop.
  #initial sample to default to when rejecting

  prop_particle_index = sample([1:current_MCMC_state.current_sample.num_particles])
  prev_sample_params = current_MCMC_state.current_sample.params[:,prop_particle_index]
  prev_sample_state = current_MCMC_state.current_sample.state[:,prop_particle_index]
  chain_len = current_MCMC_state.burnin_len + current_MCMC_state.thin_len*current_MCMC_state.num_samples_desired
  record_index = 0
  num_acc = 0

  prop_sample_params = Array(Float64, par_dim)
  prop_sample_state = Array(Int64, state_dim)
  range_num_part = [1:num_particles]
  if do_kde
    kde_kernel = Normal(0,bandwidth)
  end
  for chainstep in 1:chain_len
    #param proposal from previous stage
    prop_particle_index = sample(range_num_part)
    prop_sample_params = samples_post_prev_stage.params[:,prop_particle_index]
    if(do_kde)
      prop_sample_params = prop_sample_params.*exp(rand(kde_kernel, par_dim))#do kde in log space
    end
    prop_sample_state = samples_post_prev_stage.state[:,prop_particle_index]

    #add perturbation as if sampling from a KDE
    prop_sample_params = prop_sample_params

    #state proposal from fwd simulation
    inside_sampler = true
    prop_sample_state =

    #accept if A > 1 or A > unif, i.e. log>0 or log>log(unif)
    log_acc_rat =
      log_measurement_density(prop_sample_state, noise_distribution, obs_molecule_index, d_obs_this_stage) -
      log_measurement_density(prev_sample_state, noise_distribution, obs_molecule_index, d_obs_this_stage)
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

  push!(current_MCMC_state.num_acc, num_acc)
  current_MCMC_state.stage = current_MCMC_state.stage + 1
  return current_MCMC_state
end

function plot_from_MCMC(current_MCMC_state, ground_truth=[])
  num_particles = current_MCMC_state.current_sample.num_particles
  par_dim = current_MCMC_state.current_sample.par_dim
  state_dim = current_MCMC_state.current_sample.state_dim

  skip_len = ceil(Int64, num_particles/1000)
  plot_from_indices = 1:skip_len:num_particles

  posterior_plots = table(1,par_dim+state_dim)
  for i in 1:(par_dim + state_dim)
    if i > par_dim
      to_plot = current_MCMC_state.current_sample.state[i,plot_from_indices]
    else
      to_plot = current_MCMC_state.current_sample.params[i,plot_from_indices]
    end
    posterior_plots[1,i] = FramedPlot(title=string("Marginal", i, "of MCMC output"))
    add(posterior_plots[1,i], Histogram(hist(vec(to_plot))...))
    if !isempty(ground_truth)
      #Add a line of the right height
      line_height = maximum(hist(vec(to_plot))[2])
      add(posterior_plots[1,i], Curve([ground_truth[i],ground_truth[i]], [line_height,0], "color", "red"))
    end
  end

  return posterior_plots
end

println("----------------------------------")

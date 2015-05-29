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
using Distributions
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
    new_object = new(params,
                                   state,
                                   size(params)[1], #par_dim
                                   size(state)[1],  #state_dim
                                   max(size(params)[2], size(state)[2]) #num_particles
                                   )
    Sample_state_and_params_type_data_check(new_object)
    return new_object
  end

  function Sample_state_and_params_type(num_particles::Int64, par_dim::Int64, state_dim::Int64)
    new_object = new(zeros(Float64, par_dim, num_particles),
                                   zeros(Int64, state_dim, num_particles),
                                   par_dim,
                                   state_dim,
                                   num_particles)
    Sample_state_and_params_type_data_check(new_object)
    return new_object
  end
end
    function Sample_state_and_params_type_data_check(ssp::Sample_state_and_params_type)
      mismatch1 = (ssp.par_dim != size(ssp.params)[1])
      mismatch2 = (ssp.state_dim != size(ssp.state)[1])
      mismatch3 = (ssp.num_particles != size(ssp.params)[2])
      mismatch4 = (ssp.num_particles != size(ssp.state)[2])
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
  burnin_len::Int64
  thin_len::Int64
  bandwidth::Float64
  do_kde::Bool
  num_acc::Array{Int64, 1}
  stage::Int64
  fwd_sim::Function
  emission_logden::Function
  current_sample::Sample_state_and_params_type
end
    MCMC_state() = MCMC_state(
      1000,    #burnin_len
      5,       #thin_len
      0.001,   #bandwidth
      false,   #do_kde
      Int64[], #num_acc
      0,       #stage
      (x -> x), #fwd_sim
      (x -> 0), #emission_logden
      Sample_state_and_params_type(0, 0, 0)
      )


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
function pMCMC!(d_obs, t_obs, MCS::MCMC_state)
  println("WARNING: pMCMC! modifies the MCMC_state objects that it is given.")
  I = length(d_obs)
  for stage = 1:I #by stage, I mean how much data has been conditioned upon. At stage 2, we've conditioned on 2 data points.
    #report progress
    println(string("In pMCMC at stage ", stage, "."))
    #check data to see if pMCMC_single_stage! introduced any errors
    Sample_state_and_params_type_data_check(MCS.current_sample)

    #fold in more data
    if stage>1
      T_sim_this_stage = t_obs[stage] - t_obs[stage]
    else
      T_sim_this_stage = t_obs[stage]
    end
    d_obs_this_stage = d_obs[stage]

    #Run the sampler
    pMCMC_single_stage!(d_obs_this_stage, T_sim_this_stage, MCS)
  end
  return MCS
end

#This subroutine converts a large sample from P(params, state at time i-1|data to time i-1) into
#a large sample from P(params, state at time i|data to time i)
function pMCMC_single_stage!(d_obs_this_stage, T_sim_this_stage, MCS::MCMC_state)
  #new empty array
  samples_for_next_stage = Sample_state_and_params_type(MCS.current_sample.num_particles,MCS.current_sample.par_dim, MCS.current_sample.state_dim)

  #------------------------------Loop setup------------------------------

  #Initial values
  prop_particle_index = sample([1:MCS.current_sample.num_particles])
  prev_sample_params = MCS.current_sample.params[:,prop_particle_index]
  prev_sample_state = MCS.current_sample.state[:,prop_particle_index]
  #empty arrays for proposals
  prop_sample_params = Array(Float64, MCS.current_sample.par_dim)
  prop_sample_state = Array(Int64, MCS.current_sample.state_dim)

  chain_len = MCS.burnin_len + MCS.thin_len*MCS.current_sample.num_particles #says how many iterations we'll do
  record_index = 0 #Says how many samples we've saved.
  num_acc = 0 #Says how many proposals we've accepted.
  range_num_part = [1:MCS.current_sample.num_particles] #Helps with fast resampling if I allocate this outside the loop.
  if MCS.do_kde
    kde_kernel = Normal(0,MCS.bandwidth)
  end

  #------------------------------Loop------------------------------
  for chainstep in 1:chain_len
    #param proposal from previous stage
    prop_particle_index = sample(range_num_part)
    prop_sample_params = MCS.current_sample.params[:,prop_particle_index]
    prop_sample_state = MCS.current_sample.state[:,prop_particle_index]
    if(MCS.do_kde)
      prop_sample_params = prop_sample_params.*exp(rand(kde_kernel, MCS.current_sample.par_dim))#do kde as if in log space
    end

    #This LF-MCMC requires you to generate the proposal for the hidden state by simulating, conditioned on the proposed parameters.
    prop_sample_state = MCS.fwd_sim(prop_sample_state, prop_sample_params, T_sim_this_stage)

    #accept if A > 1 or A > unif, i.e. log>0 or log>log(unif)
    log_acc_rat =
      MCS.emission_logden(prop_sample_state, d_obs_this_stage) -
      MCS.emission_logden(prev_sample_state, d_obs_this_stage)
    if (log_acc_rat > 0) || (log_acc_rat > log(rand(1)[1]))
      acc_sample_state = prop_sample_state
      acc_sample_params = prop_sample_params
      num_acc = num_acc + 1
    else
      acc_sample_state = prev_sample_state
      acc_sample_params = prev_sample_params
    end

    #record after burnin unless thinned out
    if ((chainstep>MCS.burnin_len) && (chainstep%MCS.thin_len==1))
      record_index = record_index + 1
      samples_for_next_stage.params[:,record_index] = acc_sample_params
      samples_for_next_stage.state[:,record_index] = acc_sample_state
    end

    #Make sure the next round compares to the right items
    prev_sample_state = acc_sample_state
    prev_sample_params = acc_sample_params
  end

  #In place of return statement, this function modifies its args.
  push!(MCS.num_acc, num_acc)
  MCS.stage = MCS.stage + 1
  MCS.current_sample = samples_for_next_stage
end

using Winston

function plot_save_marginals(MCS::pMCMC_julia.MCMC_state, save_folder::String, ground_truth_params=[], ground_truth_state=[])
  num_particles = MCS.current_sample.num_particles
  par_dim = MCS.current_sample.par_dim
  state_dim = MCS.current_sample.state_dim

  skip_len = ceil(Int64, num_particles/1000)
  plot_from_indices = 1:skip_len:num_particles

  for i in 1:par_dim
    to_plot = MCS.current_sample.params[i,plot_from_indices]
    posterior_plot = plot(Int64[], Int64[], title=string("Marginal", i, "of MCMC output (pars, log scale)"))
    add(posterior_plot, Histogram(hist(vec(log10(to_plot)), 10)...))
    if !isempty(ground_truth_params)
      #Add a line of the right height
      line_height = maximum(hist(vec(log10(to_plot)), 10)[2])
      add(posterior_plot, Curve([log10(ground_truth_params[i]),log10(ground_truth_params[i])], [line_height,0], "color", "red"))
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

end

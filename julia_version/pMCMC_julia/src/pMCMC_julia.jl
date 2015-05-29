# This module implements the particle MCMC method from Darren Wilkinson's 2009 paper
# "Parameter inference for stochastic kinetic models of bacterial gene regulation: a Bayesian approach to systems biology"
# The method is for Bayesian inference on hidden markov models of continuous time-series data.
# This code currently supports only state spaces consisting of ordered lists of integers and parameter spaces equal to R^n.
# The main call is like this:


module pMCMC_julia
using Distributions

#---------------------------------MCMC_in_out_prefs--------------------------------------------
#Holds samples, inputs and preferences for the pMCMC algo
type MCMC_state
  #in
    t_obs::Array{Float64, 1}
    d_obs::Array{Int64, 1}

  #out/state
    num_acc::Array{Int64, 1}
    stage::Int64
    param_sample::Array{Float64,2}
    state_sample::Array{Int64,2}
    par_dim::Int64
    state_dim::Int64
    num_particles::Int64

  #prefs
    burnin_len::Int64
    thin_len::Int64
    bandwidth::Float64
    do_kde::Bool
    fwd_sim::Function
    emission_logden::Function
end

MCMC_state() = MCMC_state(
  #in
    Float64[], #t_obs
    Int64[],   #d_obs

  #out/state
    Int64[], #     num_acc
    0,       #     stage
    zeros(Float64, 0, 0),  # param_sample
    zeros(Int64, 0, 0),      # state_sample
    0,       #     par_dim::Int64
    0,       #     state_dim::Int64
    0,       #     num_particles::Int64

  #prefs
    1000,    #    burnin_len
    5,       #    thin_len
    0.001,   #    bandwidth
    false,   #    do_kde
    (x -> x), #    fwd_sim
    (x -> 0), #    emission_logden
)

function MCMC_state_print(MCS)
  println("t_obs:", MCS.t_obs)
  println("d_obs:", MCS.d_obs)
  println("num_acc:", MCS.num_acc)
  println("size(param_sample):", size(MCS.param_sample))
  println("size(state_sample):", size(MCS.state_sample))
  println("par_dim:", MCS.par_dim)
  println("state_dim:", MCS.state_dim)
  println("num_particles:", MCS.num_particles)
  println("burnin_len:", MCS.burnin_len)
  println("thin_len:", MCS.thin_len)
  println("bandwidth:", MCS.bandwidth)
  println("do_kde:", MCS.do_kde)
  println("fwd_sim:", MCS.fwd_sim)
  println("emission_logden:", MCS.emission_logden)
end

function MCMC_state_data_check(MCS::MCMC_state)
  mismatch = Bool[];
  push!(mismatch, length(MCS.t_obs)!= length(MCS.d_obs))
  push!(mismatch, size(MCS.param_sample)[1] == MCS.par_dim)
  push!(mismatch, size(MCS.param_sample)[2] == MCS.num_particles)
  push!(mismatch, size(MCS.state_sample)[1] == MCS.state_dim)
  push!(mismatch, size(MCS.state_sample)[2] == MCS.num_particles)
  for(problem in mismatch)
    if(problem)
      println("mismatch in MCMC_state object. Printing.")
      MCMC_state_print(MCS)
      error("mismatch in MCMC_state object. See printout.")
    end
  end
end


#This function loops over the dataset. At iteration i, it calls an MCMC-based
#subroutine to convert a large sample from P(params|data to time i-1) into
#a large sample from P(params|data to time i)
function pMCMC!(MCS::MCMC_state)
  println("WARNING: pMCMC! modifies the MCMC_state objects that it is given.")
  I = length(MCS.d_obs)
  for stage = 1:I #by stage, I mean how much data has been conditioned upon. At stage 2, we've conditioned on 2 data points.
    #report progress
    println(string("In pMCMC at stage ", stage, "."))

    #check data to see if pMCMC_single_stage! introduced any errors
    MCMC_state_data_check(MCS)

    #save the current state
    println("Still need to finish code to save intermediate data in pMCMC!().")

    #fold in more data
    if stage>1
      T_sim_this_stage = MCS.t_obs[stage] - MCS.t_obs[stage]
    else
      T_sim_this_stage = MCS.t_obs[stage]
    end
    d_obs_this_stage = MCS.d_obs[stage]

    #Run the sampler
    pMCMC_single_stage!(d_obs_this_stage, T_sim_this_stage, MCS)
  end
  return MCS
end

#This subroutine converts a large sample from P(params, state at time i-1|data to time i-1) into
#a large sample from P(params, state at time i|data to time i)
function pMCMC_single_stage!(d_obs_this_stage, T_sim_this_stage, MCS::MCMC_state)
  #------------------------------Loop setup------------------------------

  #Initial values
  prop_particle_index = sample([1:MCS.num_particles])
  prev_sample_params = MCS.params[:,prop_particle_index]
  prev_sample_state = MCS.state[:,prop_particle_index]

  #empty arrays for individual proposals
  prop_sample_params = zeros(Float64, MCS.par_dim)
  prop_sample_state = zeros(Int64, MCS.state_dim)

  #empty arrays for all accepted samples from this stage
  state_sample_for_next_stage = zeros(Float64, MCS.par_dim, MCS.num_particles)
  param_sample_for_next_stage = zeros(Float64, MCS.state_dim, MCS.num_particles)

  chain_len = MCS.burnin_len + MCS.thin_len*MCS.num_particles #says how many iterations we'll do
  record_index = 0 #Says how many samples we've saved.
  num_acc = 0 #Says how many proposals we've accepted.
  range_num_part = [1:MCS.num_particles] #Helps with fast resampling if I allocate this outside the loop.
  if MCS.do_kde
    kde_kernel = Normal(0,MCS.bandwidth)
  end

  #------------------------------Loop------------------------------
  for chainstep in 1:chain_len
    #param proposal from previous stage
    prop_particle_index = sample(range_num_part)
    prop_sample_params = MCS.params[:,prop_particle_index]
    prop_sample_state = MCS.state[:,prop_particle_index]
    if(MCS.do_kde)
      prop_sample_params = prop_sample_params.*exp(rand(kde_kernel, MCS.par_dim))#do kde as if in log space
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
  MCS.state_sample = state_sample_for_next_stage
  MCS.param_sample = param_sample_for_next_stage
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

module pMCMC_julia
using Distributions
using HDF5, JLD

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

#Checks an object of type Sample_state_and_params_type_data_check for internal consistency.
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
  verbose::Bool
  pause_len::Int64
  burnin_len::Int64
  thin_len::Int64
  bandwidth_multiplier::Float64
  bw_max::Float64
  bw_min::Float64
  do_kde::Bool
  num_acc::Array{Int64, 1}
  stage::Int64
  fwd_sim::Function
  emission_logden::Function
  current_sample::Sample_state_and_params_type
  save_intermed::Bool
  save_path::String
end

###################
#constructor to avoid breaking backwards compatibility
MCMC_state(
  verbose::Bool,
  pause_len::Int64,
  burnin_len::Int64,
  thin_len::Int64,
  bandwidth_multiplier::Float64,
  bw_max::Float64,
  bw_min::Float64,
  do_kde::Bool,
  num_acc::Array{Int64, 1},
  stage::Int64,
  fwd_sim::Function,
  emission_logden::Function,
  current_sample::Sample_state_and_params_type,
  #save intermed is lacking here but not below
  save_path::String) = MCMC_state(

  verbose::Bool,
  pause_len::Int64,
  burnin_len::Int64,
  thin_len::Int64,
  bandwidth_multiplier::Float64,
  bw_max::Float64,
  bw_min::Float64,
  do_kde::Bool,
  num_acc::Array{Int64, 1},
  stage::Int64,
  fwd_sim::Function,
  emission_logden::Function,
  current_sample::Sample_state_and_params_type,
  true::Bool,
  save_path::String)
###################

#empty default constructor
MCMC_state() = MCMC_state(
  false,   #verbose
  3,       #pause_len
  1000,    #burnin_len
  5,       #thin_len
  1,       #bandwidth_multiplier
  1,       #bw_max
  0.1,     #bw_min
  true,    #do_kde
  Int64[], #num_acc
  0,       #stage
  (x -> x), #fwd_sim
  (x -> 0), #emission_logden
  Sample_state_and_params_type(0, 0, 0),
  false,     #save_intermed::Bool
  "" #save_path
)

#Saves the MCMC_state object to folder MCS.save_path, using the file name save_file.
function MCS_save(save_name::String, MCS::MCMC_state)
  temp = MCS.save_path
  save_to = string(temp,"/", save_name)
  println("MCS_save is writing to: ", save_to)
  save(save_to,
        "verbose", MCS.verbose,
        "pause_len", MCS.pause_len,
        "burnin_len", MCS.burnin_len,
        "thin_len", MCS.thin_len,
        "bandwidth_multiplier", MCS.bandwidth_multiplier,
        "bw_max", MCS.bw_max,
        "bw_min", MCS.bw_min,
        "do_kde", MCS.do_kde,
        "num_acc", MCS.num_acc,
        "stage", MCS.stage,
        "current_sample_params", MCS.current_sample.params,
        "current_sample_state", MCS.current_sample.state,
        "save_intermed",MCS.save_intermed,
        "save_path", MCS.save_path)
end

#Loads the MCMC_state object from the file at save_path
function MCS_load(save_path)
   loaded_dict = load(save_path)
  MCS = MCMC_state(loaded_dict["verbose"],
                  loaded_dict["pause_len"],
                  loaded_dict["burnin_len"],
                  loaded_dict["thin_len"],
                  loaded_dict["bandwidth_multiplier"],
                  loaded_dict["bw_max"],
                  loaded_dict["bw_min"],
                  loaded_dict["do_kde"],
                  loaded_dict["num_acc"],
                  loaded_dict["stage"],
                  x->x, #fwd_sim
                  x->0, #emission_logden
                  Sample_state_and_params_type(loaded_dict["current_sample_params"],loaded_dict["current_sample_state"]),
                  loaded_dict["save_intermed"],
                  loaded_dict["save_path"])
  return MCS
end

#Checks for internal inconsistency in an MCMC_state object.
#Not thorough objects with inconsistent data may escape detection.
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


    #save samples
    if MCS.save_intermed
      println("In pMCMC!, saving to ", MCS.save_path)
      pMCMC_julia.MCS_save(string("stage",MCS.stage, "sample"), MCS)
    end

    #fold in more data
    if stage>1
      T_sim_this_stage = t_obs[stage] - t_obs[stage]
    else
      T_sim_this_stage = t_obs[stage]
    end
    d_obs_this_stage = d_obs[stage]
    if MCS.verbose
      println("T_sim_this_stage is: ", T_sim_this_stage)
      println("d_obs is currently: ", d_obs_this_stage)
    end
    #Run the sampler
    pMCMC_single_stage!(d_obs_this_stage, T_sim_this_stage, MCS)
  end
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
  if MCS.verbose
    println("prop_particle_index is: ", prop_particle_index)
    println("prev_sample_params is: ", prev_sample_params)
    println("prev_sample_state is: ", prev_sample_state)
  end

  #empty arrays for proposals
  prop_sample_params = Array(Float64, MCS.current_sample.par_dim)
  prop_sample_state = Array(Int64, MCS.current_sample.state_dim)

  chain_len = MCS.burnin_len + MCS.thin_len*MCS.current_sample.num_particles #says how many iterations we'll do
  record_index = 0                                                           #Says how many samples we've saved.
  num_acc = 0                                                                #Says how many proposals we've accepted.
  range_num_part = [1:MCS.current_sample.num_particles] #Helps with fast resampling if I allocate this outside the loop.

  #set up the KDE
  silverman_bw = ones(Float64, MCS.current_sample.par_dim)
  for d in 1:MCS.current_sample.par_dim
    silverman_bw[d] = MCS.current_sample.num_particles^(-1/(4+MCS.current_sample.par_dim))*std(log(MCS.current_sample.params[d, :]))
  end
  bw_this_stage = MCS.bandwidth_multiplier*silverman_bw
  bw_this_stage = [minimum([cand_bw, MCS.bw_max]) for cand_bw in bw_this_stage]
  bw_this_stage = [maximum([cand_bw, MCS.bw_min]) for cand_bw in bw_this_stage]
  bw_this_stage = mean(bw_this_stage)
  #println("Current bandwidths: ", bw_this_stage)
  if MCS.do_kde
    kde_kernel = Normal(0,1)
  end
  if MCS.verbose
    println("bw_this_stage is: ", bw_this_stage)
    println("kde_kernel is: ", kde_kernel)
  end
  #------------------------------Loop------------------------------
  for chainstep in 1:chain_len
    #param proposal from previous stage
    prop_particle_index = sample(range_num_part)
    prop_sample_params = MCS.current_sample.params[:,prop_particle_index]
    prop_sample_state = MCS.current_sample.state[:,prop_particle_index]
    if MCS.verbose
      println("prop_particle_index is: ", prop_particle_index)
      println("prop_sample_params before KDE is: ", prop_sample_params)
      println("prop_sample_state before fwd sim is: ", prop_sample_state)
    end

    #perturb samples additive in log-space according to KDE kernel
    if(MCS.do_kde)
      prop_sample_params = prop_sample_params.*exp(bw_this_stage.*rand(kde_kernel, MCS.current_sample.par_dim))#do kde as if in log space
    end
    #This LF-MCMC requires you to generate the proposal for the hidden state by simulating, conditioned on the proposed parameters.
    prop_sample_state = MCS.fwd_sim(prop_sample_state, prop_sample_params, T_sim_this_stage, MCS.verbose)
    if MCS.verbose
      println("prop_sample_params after KDE is: ", prop_sample_params)
      println("prop_sample_state after fwd sim is: ", prop_sample_state)
    end


    #Compute acceptance ratio and unif for comparison
    logunif = log(rand(1)[1])
    log_acc_rat =
      MCS.emission_logden(prop_sample_state, d_obs_this_stage) -
      MCS.emission_logden(prev_sample_state, d_obs_this_stage)
    if MCS.verbose
      println("logden(prop_sample_state, d_obs_this_stage) is: ", MCS.emission_logden(prop_sample_state, d_obs_this_stage))
      println("logden(prev_sample_state, d_obs_this_stage) is: ", MCS.emission_logden(prev_sample_state, d_obs_this_stage))
      println("log_acc_rat is: ", log_acc_rat)
      println("logunif is: ", logunif)
    end

    #accept if A > 1 or A > unif, i.e. logA>0 or logA>log(unif)
    if (log_acc_rat > 0) || (log_acc_rat > logunif)
      acc_sample_state = prop_sample_state
      acc_sample_params = prop_sample_params
      num_acc = num_acc + 1
      if MCS.verbose
        println("Accepted!")
      end
    else
      acc_sample_state = prev_sample_state
      acc_sample_params = prev_sample_params
      if MCS.verbose
        println("Rejected!")
      end
    end
    if MCS.verbose
      println("acc_sample_params is: ", acc_sample_params)
      println("acc_sample_state is: ", acc_sample_state)
      pMCMC_julia.eric_pause(MCS.pause_len)
    end

    #record after burnin unless thinned out
    if ((chainstep>MCS.burnin_len) && (chainstep%MCS.thin_len==1))
      record_index = record_index + 1
      samples_for_next_stage.params[:,record_index] = acc_sample_params
      samples_for_next_stage.state[:,record_index] = acc_sample_state
      if MCS.verbose
        println("chainstep is: ", chainstep, ". Recording as samples:" )
        println("acc_sample_params as recorded is: ", acc_sample_params)
        println("acc_sample_state as recorded  is: ", acc_sample_state)
        pMCMC_julia.eric_pause(MCS.pause_len)
      end
    end

    #Make sure the next round compares to the right items
    prev_sample_state = acc_sample_state
    prev_sample_params = acc_sample_params
  end
  if MCS.verbose
    println("prev_sample_state updated to: ", prev_sample_state)
    println("prev_sample_params updated to: ", prev_sample_params)
    pMCMC_julia.eric_pause(MCS.pause_len)
  end
  #In place of return statement, this function modifies its args.
  push!(MCS.num_acc, num_acc)
  MCS.stage = MCS.stage + 1
  MCS.current_sample = samples_for_next_stage
  if MCS.verbose
    println("MCS.stage updated to: ", MCS.stage)
    println("MCS.current_sample.params[:, 1:2] updated to: ", MCS.current_sample.params[:, 1:2])
    println("MCS.current_sample.state[:, 1:2] updated to: ", MCS.current_sample.state[:, 1:2])
    pMCMC_julia.eric_pause(MCS.pause_len)
  end
end


#Pauses for something within an order of magnitude of num_secs seconds.
#May fail for small values.
function eric_pause(num_secs)
  println("-------------pausing on the order of ", num_secs, " seconds---------------")
  for i in 1:round(num_secs*2.5e9)
    j = i^2
  end
end


#Implementation of a Geweke test for this pMCMC module using a
#standard normal noise model and a simple forward model x(t+1) = x(t) + round(theta).
#The simulation time is 1.
#The prior on theta is unif(0, 10).
#The ``prior'' on x(0) is a point mass at 1.
#first, helpers.
  geweke_noise_dist = Normal(0,1)
  geweke_prior_params() = [10*rand(1)]
  geweke_prior_state() = [1]
  geweke_fwd_sim(prop_sample_state, prop_sample_params, T_sim_this_stage, verbose) = prop_sample_state + convert(Int64, round(prop_sample_params[1]))
  geweke_emission_logden(x_current, d_obs) = logpdf(geweke_noise_dist, d_obs-x_current[1,1])
  geweke_get_d(sample_state) = sample_state[1] + rand(geweke_noise_dist)
  #Need to be able to perform a single MH iteration.
  #One easy way: imagine a run with just a few particles.
  #Might as well run 10 times and get some bang for your buck, hence the thin.


function single_MCMC_step(obs::Float64)
  num_particles = 1000
  params_2d = zeros(Float64, length(geweke_prior_params()), num_particles)
  state_2d = zeros(Int64, length(geweke_prior_state()), num_particles)
  for i in 1:num_particles
    params_2d[:,i] = geweke_prior_params()
    state_2d[:,i] = geweke_prior_state()
  end

  geweke_MCS = MCMC_state(
   false,#verbose::Bool
    1, #pause_len::Int64
    0, #burnin_len::Int64
    10, #thin_len::Int64
    0, #bandwidth_multiplier::Float64
    0, #bw_max::Float64
    0, #bw_min::Float64
    false, #do_kde::Bool
    Int64[], #num_acc::Array{Int64, 1}
    0, #stage::Int64
    geweke_fwd_sim, #fwd_sim::Function
    geweke_emission_logden, #emission_logden::Function
    Sample_state_and_params_type(params_2d, state_2d), #current_sample::Sample_state_and_params_type
    false, #save_intermed::Bool
    "Users/EricKernfeld/Desktop" #save_path::String
   )#MCMC_state
  pMCMC_single_stage!(obs, [1.0], geweke_MCS)#pMCMC
  return_particle = sample(1:geweke_MCS.current_sample.num_particles)
  return geweke_MCS.current_sample.params[:, return_particle],geweke_MCS.current_sample.state[:,return_particle]
end

function simple_geweke_test(samp_size)

  #because the pMCMC module is written for more complicated situations, the state/params arrays must be 2d.
  forward_params = zeros(Float64, 1, samp_size)
  forward_states = zeros(Int64, 1, samp_size)
  forward_obs = zeros(Float64, samp_size)

  #Geweke compares two samples that should be the same.
  #One comes from a forward sampling process: just get synthetic data from the generative model.
  for i in 1:samp_size
    forward_params[:, i] = geweke_prior_params()
    forward_states[:, i] = geweke_fwd_sim(geweke_prior_state(), forward_params[:, i], 1, false)
    forward_obs[i] = geweke_get_d(forward_states[:,i])
  end

  #The other sample that Geweke looks at comes from a Markov chain.
  #start with empty space and a single forward draw theta, x, d
  mc_params = zeros(Float64, 1,samp_size)
  mc_states = zeros(Int64, 1,samp_size)
  mc_obs = zeros(Float64, 1,samp_size)
  mc_params[:,1] = forward_params[:,1]
  mc_states[:,1] = forward_states[:,1]
  mc_obs[1] = forward_obs[1]

  #alternate the MCMC transition with regenerating data many times
  for i in 2:samp_size
    #Resample theta using pMCMC
    mc_params[:,i], mc_states[:,i] = single_MCMC_step(mc_obs[i-1])
    #Resample d using theta, x
    mc_obs[i] = geweke_get_d(mc_states[:,i])
  end

  geweke_fwd_samp = (forward_params,forward_states,forward_obs)
  geweke_mc_samp = (mc_params,mc_states,mc_obs)
  return (geweke_fwd_samp, geweke_mc_samp)
end

end


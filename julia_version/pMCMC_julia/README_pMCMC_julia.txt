README_pMCMC_julia

This module implements the particle MCMC method from Darren Wilkinson's 2009 paper
"Parameter inference for stochastic kinetic models of bacterial gene regulation: a Bayesian approach to systems biology"
The method is for Bayesian inference on hidden markov models of continuous time-series data.
This code currently supports only state spaces consisting of ordered lists of integers, parameter spaces equaling R^n, and real-valued scalar observations.

The easiest way to interface with this code is to use the Prelim_experiments module. To use it directly, an example is given below, but there are significant restrictions on each argument. The documentation explaning these restrictions refers both to the variable names and to the first example below.

Docs for members of data type MCMC_state
  verbose: a Boolean value. If true, the sampler prints every detail of its actions.
  pause_len: how long to pause upon printing in verbose mode (in seconds, but it's homegrown and may be very inaccurate). Real scalar.
  burnin_len: the number of Metropolis-Hastings steps to take before beginning to record samples. Positive integer.
  thin_len: In the Metropolis-Hastings, one of every thin_len samples gets recorded. Positive integer.
  bandwidth_multiplier: Resampling uses a KDE with the Silverman rule bandwidth, ied by this. Recommended value: 0.1. Real scalar.
  bw_max: maximum bandwidth to use when resampling. Real scalar.
  bw_min: minimum bandwidth to use when resampling. Real scalar.
  do_kde: Boolean value indicating whether to use a KDE when re-sampling.
  num_acc: Array containing the number of samples accepted in each round. Not meant to be modified.
  stage: integer indicating what stage the sampler is at. Not meant to be modified.
  fwd_sim: Function. See below.
  emission_logden: Function. See below.
  current_sample: composite type containing samples. See below.
  save_intermed: Boolean value. If true, sampler saves MCMC_state object to save_path at each stage.
  save_path: String. See docs on save_intermed.


Docs for remaining members, phrased in terms of the example code below
  large_sample_of_your_parameters--2D array of floats
  large_sample_of_your_initial_states--2D array of integers
    These two should have matching second dimensions. For either one, indexing by [:, i] should give a single sample from your prior distribution.
  The variable marked as <a simulator of your choice, function-valued> should be function-valued. It may be stochastic.
     It should take:
       parameters, (latent) initial state, desired duration, verbose
     It should return:
       the state after a simulation of the desired duration

  The variable marked as <log probability density for your emissions, function-valued> should be a function to compute the log-density of the measurement error
    It should take:
       system_state, 1d array of integers
       my_observation, a real scalar
    It should return:
       log-density of the measurement error assuming the system was in system_state and emitted my_observation

#------------------------start example---------------------------
using pMCMC_julia
your_MCS = pMCMC_julia.MCMC_state()
#your_MCS.fwd_sim = <a simulator of your choice, function-valued>
#your_MCS.emission_logden = <log probability density for your emissions, function-valued>
#your_MCS.current_sample = pMCMC_julia.sample_state_and_params(<large_sample_of_your_parameters>, <large_sample_of_your_initial_states>)
pMCMC!(your_data, your_timepoints, your_MCS)

#your_MCS.current_sample gets modified by pMCMC!.
#These lines each access a single posterior sample.
#your_MCS.current_sample.params[:, i]
#your_MCS.current_sample.state[:, i]
#------------------------end example---------------------------


Other functions in the module meant to be used as part of the interface include:

  MCS_save(save_name::String, MCS::MCMC_state)
    Saves the MCMC_state object to folder MCS.save_path, using the file name save_file. Before using this, set MCS.save_path wherever you would like the file.

  MCS_load(save_path)
    Loads and returns the MCMC_state object from the file at save_path. Not guaranteed to recover data saved by any function other than MCS_save.

  Sample_state_and_params_type(prior_params::Array{Float64, 2}, prior_state::Array{Int64, 2})
    Packs the prior samples into a specialized data type. This object will be modified when you run the pMCMC.
    For either argument, indexing by [:, i] should yield a single sample from your prior distribution.

  simple_geweke_test(samp_size)
    Performs a Geweke test of the function pMCMC!(), generating samp_size draws by MCMC and from the generative model.
    The simplest interface to this is through the script geweke_test_script.jl.
    (see also John Geweke, 2004, JASA Vol 99 No 467: "Getting It Right: Joint Distribution Tests of Posterior Simulators")


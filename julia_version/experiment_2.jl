#This file lays out an experiment where no reactions occur.
#The forward sim can handle a rate of 0, so I put in that in, plus a
#stoichiometry matrix that gives no change at each reaction event.
#The prior is x~discrete uniform on 0 to 500, with observations d~x+binomial(40, 0.5)
#and true state x=50 molecules, so the posterior should equal the likelihood, renormalized.
#It will be a product of binomial terms.

#File layout:
#set up measurement error (two functions: one to sample and another to evaluate the density)
#then load up the essentials
#then set up the other inputs (e.g. the stoichiometry matrix)

function log_measurement_density(x_current, noise_distribution, molecule_index, d_obs)
  log_density = logpdf(noise_distribution, d_obs-x_current[molecule_index])
  return log_density
end

function do_measurement_error(x_current, noise_distribution,molecule_index)
  d_obs = x_current[molecule_index] + rand(noise_distribution, 1)[1]
  return d_obs
end

#Essentials: tells Julia where to find the forward simulation,
#the MCMC routines, the Distributions package, and a custom data type that I made.
include("prelim_software_setup.jl")

# This function sets up the following arguments for pMCMC() and the forward simulation.
# d_obs -- drawn via fwd simulation using the true params.
# t_obs -- a regular grid chosen by the experimenter
# rxn_times -- drawn via fwd simulation using the true params.
# prior_sample -- drawn from a vague prior.
# num_samples_desired -- 1e6
# burnin_len -- 1e3,
# thin_len -- 5,
# A system with one molecule type. These will cause exponential decay (left entry)
# and constant spontaneous generation (right entry).
#   sto_mat -- [-1, 1]
#   rxn_entry_mat -- [1,0]
#   true_rxn_rates -- [0.001, 0.1]
# By convention, I set things up so that the first molecule is the one measured.
# For experiments with Wilkinson's examples, have to consider this again.
#   molecule_index -- 1
function pMCMC_setup_experiment_2(;actually_run_MCMC=false)
  true_rxn_rates = [0.0]
  true_init_x = Int64[50]
  sto_mat = zeros(Int64, 1,1)
  rxn_entry_mat = zeros(Int64, 1,1)
  t_interval = 300
  num_intervals = 24
  T_sim = t_interval*num_intervals
  t_obs = t_interval*[1:num_intervals]
  molecule_index = 1
  noise_distribution = Binomial(40, 0.5)
  inside_sampler = false

  bandwidth = 1
  do_kde = false
  num_samples_desired = 100000
  burnin_len = 1e3
  thin_len = 5

  #-----------------Forward simulation to generate data------------------
  x_path, current_x, num_rxns_occ, rxn_types,rxn_times, t_spent = gillespie(true_init_x,
                                                                         sto_mat,
                                                                         rxn_entry_mat,
                                                                         true_rxn_rates,
                                                                         T_sim,
                                                                         inside_sampler)
  x_index_matching_t_obs = Array(Int64, length(t_obs))
  x_obs = Array{Int64, 1}[] #Noiseless verion of observed values
  d_obs = zeros(Float64,length(t_obs)) #Observed values
  for j in 1:length(t_obs)
    #Assume at first that t_obs[j] is after all rxns have happened.
    push!(x_obs, x_path[length(rxn_times)])
    #Test other possible windows between reaction times. If
    #one of them works out, exit the loop knowing that t_obs[j] includes only first i-1 rxns
    for i in 2:length(rxn_times)
      if(rxn_times[i-1]<=t_obs[j] && t_obs[j] <= rxn_times[i])
        x_index_matching_t_obs[j] = i - 1
        x_obs[j] = x_path[x_index_matching_t_obs[j]]
        break
      end
    end
    d_obs[j] = do_measurement_error(x_obs[j], noise_distribution, molecule_index)
  end

  #-----------------Inference-------------------------------
  if(actually_run_MCMC==true)
    #clairvoyant prior on reaction rate: this way the simulation won't ever change the state.
    prior_sample_thetas = zeros(2,num_samples_desired)
    prior_sample = Sample_state_and_params_type(prior_sample_thetas, reshape(sample(1:500, num_samples_desired),1,num_samples_desired))

    #Some records for later
    metadata_to_save = string("This test was run at time ", now(), " with ",
                            "observations at intervals of ", t_interval,
                            " from time zero to ", T_sim,
                            " with a noise sd of ", std(noise_distribution), ". ",
                            " There were ",  num_samples_desired, " particles, ",
                            " with a clairvoyant spike-at-zero prior on reaction rate, ",
                            " a discrete unif(0,500) prior on the state, ",
                            " and true state of ", true_init_x, ".")

    posterior_sample, num_acc = pMCMC(d_obs,
                             t_obs,
                             prior_sample,
                             num_samples_desired,
                             burnin_len,
                             thin_len,
                             noise_distribution,
                             molecule_index,
                             sto_mat,
                             rxn_entry_mat,
                             bandwidth,
                             do_kde)
    return posterior_sample, prior_sample, true_rxn_rates, t_obs, d_obs, x_obs, x_path, true_init_x, num_acc, metadata_to_save
  end
  return true_rxn_rates, t_obs, d_obs, x_obs, x_path, true_init_x
end
include("prelim_software_setup.jl")

#Set things up, bt don't pull the trigger on the pMCMC
#true_rxn_rates, t_obs, d_obs, x_obs, x_path, true_init_x = pMCMC_setup_experiment_2(actually_run_MCMC=false)

tic()
#Actually run the sampler; update metadata_to_save with time taken
  Profile.clear()
  posterior_sample,prior_sample, true_rxn_rates, t_obs, d_obs, x_obs, x_path, true_init_x, num_acc, metadata_to_save = @profile pMCMC_setup_experiment_2(actually_run_MCMC=true)
  ProfileView.view()
  time_taken = toc()
  metadata_to_save = string(metadata_to_save, " It took ", time_taken, " seconds.")
  num_acc

#make a plot of the simulated trajectory
x_obs_vals = Array(Int64,length(t_obs))
  for i in 1:length(x_obs)
    x_obs_vals[i] = x_obs[i][1]
  end
  sim_trajec_plot = FramedPlot(
           title="Hidden (black) and observed (red) states from simulation",
           xlabel="Time",
           ylabel="State")
  add(sim_trajec_plot, Points(t_obs, x_obs_vals))
  add(sim_trajec_plot, Points(t_obs, d_obs, "color","red"))
#plot prior and posterior
post_plot = FramedPlot(title="Posterior histogram")
  add(post_plot, Histogram(hist(vec(posterior_sample.state[1, :]))...))
  line_height = maximum(hist(vec(posterior_sample.state[1,:]))[2])
  add(post_plot, Curve([true_init_x[1],true_init_x[1]], [line_height,0], "color", "red"))
pri_plot = FramedPlot(title="Prior histogram")
  add(pri_plot, Histogram(hist(vec(prior_sample.state[1,:]))...))
  line_height = maximum(hist(vec(prior_sample.state[1,:]))[2])
  add(pri_plot, Curve([true_init_x[1],true_init_x[1]], [line_height,0], "color", "red"))


#Make a new folder named by experiment and date. Save profiler data, plots of prior and posterior, and

today_filepath = string("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/experiment2",now())
mkdir(today_filepath)
#to save prior and posterior plots and plot of data
  savefig(post_plot, string(today_filepath, "/post.png"))
  savefig(pri_plot, string(today_filepath, "/prior.png"))
  savefig(sim_trajec_plot, string(today_filepath, "/sim_trajec.pdf"))

#To save important metadata
  @save string(today_filepath, "/metadata") metadata_to_save

#To save profiler data
bt, lidict = Profile.retrieve()
  @save string(today_filepath, "/profiler_data") bt lidict


#To load profiler data
  #using HDF5, JLD, ProfileView
  #@load today_filepath
  #ProfileView.view(bt, lidict=lidict)


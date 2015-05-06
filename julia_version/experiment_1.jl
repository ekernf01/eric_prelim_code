#File layout:
#set up measurement error (two functions: one to sample and another to evaluate the density)
#then load up the essentials
#then set up the other inputs (e.g. the stoichiometry matrix)

function log_measurement_density(x_current, noise_sd, molecule_index, d_obs)
  log_density = logpdf(Normal(x_current[molecule_index], noise_sd), d_obs)
  return log_density
end

function do_measurement_error(x_current, noise_sd,molecule_index)
  d_obs = x_current[molecule_index] + rand(Normal(0, noise_sd), 1)[1]
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
function pMCMC_setup_experiment_1(;actually_run_MCMC=false)
  true_rxn_rates = [0.001, 0.1]
  true_init_x = [50]
  sto_mat = reshape([-1, 1],1,2)
  rxn_entry_mat = reshape([1, 0],1,2)
  t_interval = 3 #300
  num_intervals = 6 #24
  T_sim = t_interval*num_intervals
  t_obs = t_interval*[1:num_intervals]
  molecule_index = 1
  noise_sd = 100
  inside_sampler = false

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
    d_obs[j] = do_measurement_error(x_obs[j], noise_sd, molecule_index)
  end

  #-----------------Inference-------------------------------
  if(actually_run_MCMC==true)
    num_samples_desired = 10000
    burnin_len = 1e3
    thin_len = 5
    prior_sample_thetas = 10.^(-4*rand(num_samples_desired, 2)) #log uniform prior on 1*10^-4 to 1*10^0
    prior_sample = Sample_state_and_params_type(prior_sample_thetas,
                        reshape(repmat(true_init_x, num_samples_desired),num_samples_desired,1))
    posterior_sample, num_acc = pMCMC(d_obs,
                             t_obs,
                             prior_sample,
                             num_samples_desired,
                             burnin_len,
                             thin_len,
                             noise_sd,
                             molecule_index,
                             sto_mat,
                             rxn_entry_mat)
    return posterior_sample, prior_sample, true_rxn_rates, t_obs, d_obs, x_path, true_init_x, num_acc
  end
  return true_rxn_rates, t_obs, d_obs, x_obs, x_path, true_init_x
end
include("prelim_software_setup.jl")

#Set things up, bt don't pull the trigger on the takes-fucking-forever pMCMC
true_rxn_rates, t_obs, d_obs, x_obs, x_path, true_init_x = pMCMC_setup_experiment_1(actually_run_MCMC=false)

x_obs_vals = Array(Int64,length(t_obs))
for i in 1:length(x_obs)
  x_obs_vals[i] = x_obs[i][1]
end

#make a plot of the simulated trajectory
sim_trajec_plot = FramedPlot(
           title="Hidden (black) and observed (red) states from simulation",
           xlabel="Time",
           ylabel="State")
  add(sim_trajec_plot, Points(t_obs, x_obs_vals))
  add(sim_trajec_plot, Points(t_obs, d_obs, "color","red"))

#Actually run the sampler
posterior_sample,prior_sample, true_rxn_rates, t_obs, d_obs, x_path, true_init_x, num_acc = pMCMC_setup_experiment_1(actually_run_MCMC=true)
num_acc

post_plot = plot_from_MCMC(posterior_sample, 24)
pri_plot = plot_from_MCMC(prior_sample, 0)
posterior_sample.state[1:10,1]
posterior_sample.params[1:10,:]

state_hist = FramedPlot(
  title="Sampled num molecules at end of simulation",
         xlabel="Num molecules",
         ylabel="Y axis")
  add(state_hist, Histogram(hist(convert(Array{FloatingPoint,1}, posterior_sample.state[:,1]))...))


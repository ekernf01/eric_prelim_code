#File layout:
#set up measurement error (two functions: one to sample and another to evaluate the density)
#then load up the essentials
#then set up the other inputs (e.g. the stoichiometry matrix)

function log_measurement_density(x_current, noise_sd, molecule_index, d_obs)
  log_density = logpdf(Normal(x_current[molecule_index], noise_sd), d_obs)
  return log_density
end

function do_measurement_error(x_current, noise_sd)
  d_obs = x_current + rand(Normal(0, noise_sd), 1)
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
  true_init_x = 50
  sto_mat = [-1, 1]
  rxn_entry_mat = [1,0]
  T_sim = 300*24
  molecule_index = 1
  noise_sd = 0.5 #Let's start off easy

  #-----------------Forward simulation to generate data------------------
  x_path, current_x, num_rxns_occ, rxn_types,rxn_times, t_spent = gillespie(true_init_x,
                                                                         sto_mat,
                                                                         rxn_entry_mat,
                                                                         true_rxn_rates,
                                                                         T_sim)
  for j in 1:length(t_obs)
    for i in 2:length(rxn_times)
      if(rxn_times[i-1]<=t_obs[j] && t_obs[j] <= rxn_times[i])
        x_index_matching_t_obs[j] = i - 1
        break
      end
    end
    d_obs[j] = do_measurement_error(x_current=x_path[x_index_matching_t_obs[j]], noise_sd)
  end

  #-----------------Inference-------------------------------
  if(actually_run_MCMC==true)
    num_samples_desired = 1e3
    burnin_len = 1e3
    thin_len = 5
    prior_sample_thetas = 10^(-4*rand(num_samples_desired, 2)) #log uniform prior on 1*10^-4 to 1*10^0
    prior_sample = Sample_state_and_params_type(prior_sample_thetas,
                                                repmat(true_init_x, num_samples_desired))
    posterior_sample = pMCMC(d_obs,
                             t_obs,
                             prior_sample,
                             num_samples_desired,
                             burnin_len,
                             thin_len,
                             obs_molecule_index,
                             sto_mat,
                             rxn_entry_mat)
  end
end


pMCMC_setup_experiment_1(actually_run_MCMC=false)

#If you ever wonder whether doing
#
# samples_post_prev_stage = prior_sample_thetas
#
# copies by reference, well, it doesn't: run the following demo.


# num_samples_desired = 10
# true_init_x = [50]
# prior_sample_thetas = 10.^(-4*rand(num_samples_desired, 2)) #log uniform prior on 1*10^-4 to 1*10^0
#     prior_sample = Sample_state_and_params_type(prior_sample_thetas,
#                         reshape(repmat(true_init_x, num_samples_desired),num_samples_desired,1))

# samples_post_prev_stage = prior_sample_thetas
# par_dim = 2
# state_dim = 1
# samples_post_current_stage = Sample_state_and_params_type(num_samples_desired,par_dim, state_dim)

# samples_post_current_stage
# prior_sample
#Suppose you start off a test, but tell it to stop early. This will let you pick it up again.
#This function will start at the last stage you did, and it will do
#num_stages_to_try more of them.
#pMCMC has no native ability to do only the first few stages. You have to trick it:
#feed it only a segment of the data and feed it previous posterior samples as prior samples.
function pick_up_test(save_path, num_stages_to_try)
  MCS = pMCMC_julia.MCS_load(joinpath(save_path, "after_inf_sample"))

  metadata_path = joinpath(save_path, "everything_just_before_inference")
  loaded_dict = load(metadata_path)
    ep = loaded_dict["ep"]
    wilk_cri = loaded_dict["wilk_cri"]
    sim_results = loaded_dict["sim_results"]

  noise_distribution = Normal(0, ep.noise_sd)
  T_sim = ep.t_interval*ep.num_intervals
  t_obs = ep.t_interval*[1:ep.num_intervals]

  #Have to set up the model for emissions again because I can't figure out how to save functions
  function emission_logden(x_current, d_obs)
    log_density = logpdf(noise_distribution, d_obs-x_current[ep.obs_mol_ind])
    return log_density
  end
  MCS.emission_logden = emission_logden

  function fwd_sim(prev_sample_state::Array{Int64, 1}, prop_sample_params::Array{Float64, 1}, T_sim::Float64)
    for i in 1:ep.num_unks
      wilk_cri.rxn_rates[ep.unk_inds[i]] = prop_sample_params[i]
    end
    wilk_cri.init_amts = prev_sample_state
    return Chem_rxn_tools.gillespie(wilk_cri, T_sim, true)
  end
  MCS.fwd_sim = fwd_sim

#-----------------------------do the inference-------------------------------
tic()
  #You start inference. MCS.stage==0. You condition on data point #1 and run through the sampler once.
  #Then MCS.stage gets set to 1. Repeat. Thus, if MCS.stage==i, you need to add datum i+1 next.
  data_range = (MCS.stage+1):(MCS.stage+num_stages_to_try)
  println("In pick_up_test, data_range is: ", data_range)
  println("In pick_up_test, d_obs is: ", sim_results.d_obs)
  pMCMC_julia.pMCMC!(sim_results.d_obs[data_range], sim_results.t_obs[data_range], MCS)
  ep.time_taken = toc()

#-----------------------------plot/save/print the results------------------------------
  pMCMC_julia.MCS_save("after_inf_sample", MCS)

  Prelim_experiments.print_exp_prefs(ep)
  println(string("By round, the acceptance rate was ", MCS.num_acc/(MCS.burnin_len + MCS.thin_len*MCS.current_sample.num_particles)))
  make_all_plots(MCS.save_path)
end

module Prelim_experiments
  using Chem_rxn_tools
  using Distributions
  using pMCMC_julia
  using Dates
  using HDF5, JLD

type Exp_prefs
  num_samples_desired::Int64
  bw_max::Float64
  bw_min::Float64
  bandwidth_multiplier::Float64
  noise_sd::Float64
  cri_source::String
  obs_mol_name::String
  unk_names::Array{String,1}
  t_interval::Float64
  num_intervals::Array{Float64,1}
  save_folder::String

  time_of_test::DateTime
  time_taken::Float64
end

function Exp_prefs(num_samples_desired,
  bw_max,
  bw_min,
  bandwidth_multiplier,
  noise_sd,
  cri_source,
  obs_mol_name,
  unk_names,
  t_interval,
  num_intervals,
  save_folder)
  return (Exp_prefs(num_samples_desired,
          bw_max,
          bw_min,
          bandwidth_multiplier,
          noise_sd,
          cri_source,
          obs_mol_name,
          unk_names,
          t_interval,
          num_intervals,
          save_folder,
          now(),
          0))
end

function print_exp_prefs(ep::Exp_prefs)
  println("num_samples_desired: ", ep.num_samples_desired)
  println("bw_max: ", ep.bw_max)
  println("bw_min: ", ep.bw_min)
  println("bandwidth_multiplier: ", ep.bandwidth_multiplier)
  println("cri_source: ", ep.cri_source)
  println("obs_mol_name: ", ep.obs_mol_name)
  println("unk_names: ", ep.unk_names)
  println("t_interval: ", ep.t_interval)
  println("num_intervals: ", ep.num_intervals)
  println("save_folder: ", ep.save_folder)
  println("time_of_test: ", ep.time_of_test)
  println("time_taken: ", ep.time_taken)
end

function run_test_generic(ep::Exp_prefs)

  noise_distribution = Normal(0, noise_sd)
  T_sim = t_interval*num_intervals
  t_obs = t_interval*[1:num_intervals]

  #-----------------------------Load in info from SBML shorthand, set priors------------------------------
  if cri_source == "wilk"
    SBML_file = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/wilkinson_rxns_SBML_shorthand.txt"
    wilk_cri = Chem_rxn_tools.SBML_read(SBML_file)
    mols_to_show = ["SigD", "Hag", "CodY"]
    #log uniform prior on 1*10^-4 to 1*10^0 for the 3 unknowns:
      #log(kSigDprod) ∼ Unif(log{0.01}, log{100}),
      #log(kflacherep) ∼ Unif(log{0.0002}, log{2}),
      #log(kflacheunrep) ∼ Unif(log{0.001}, log{10}).
      param_sample = zeros(3,num_samples_desired)
      param_sample[1, :] = 10.^(2 - 4*rand(1,num_samples_desired))
      param_sample[2, :] = 2*(10.^(-4*rand(1,num_samples_desired)))
      param_sample[3, :] = 10.^(1 - 4*rand(1,num_samples_desired))
      state_sample = repmat(wilk_cri.init_amts,1,num_samples_desired)
  else
    wilk_cri = Chem_rxn_tools.make_demo_cri_v2()
  end

#-----------------------------Extract info from cri about unknowns------------------------------
  unk_inds = Int64[]
  num_unks = length(unk_names)
  for word in unk_names
    push!(unk_inds, Chem_rxn_tools.get_rate_indices(wilk_cri, word))
    push!(unk_rates, wilk_cri.rxn_rates[unk_inds[end]])
  end
  #Sanity check:
  for i in 1:3
    println(string(unk_names[i], " with rate ", unk_rates[i], " should be: ", wilk_cri.rxn_labels[unk_inds[i]], " with rate ", wilk_cri.rxn_rates[unk_inds[i]] ,"."))
  end

#-----------------------------Set up the sampler prefs-------------------------------

  MCS = pMCMC_julia.MCMC_state()

#BW gets set to 1 evey time
  MCS.bandwidth_multiplier = bandwidth_multiplier
  MCS.bw_max = bw_max
  MCS.bw_min = bw_min
  MCS.do_kde = true

  MCS.burnin_len = 1e3
  MCS.thin_len = 5
  MCS.current_sample = pMCMC_julia.Sample_state_and_params_type(param_sample, state_sample)

  obs_mol_ind = Chem_rxn_tools.get_chem_indices(wilk_cri, obs_mol_name)
  #Set up the model for emissions
  function emission_logden(x_current, d_obs)
    log_density = logpdf(noise_distribution, d_obs-x_current[obs_mol_ind])
    return log_density
  end
  MCS.emission_logden = emission_logden

  println("WARNING: fwd_sim WILL MODIFY wilk_cri WHEN pMCMC_julia.pMCMC IS CALLED.")
  println("It speeds things up not to copy the whole object.")
  println("Likely consequence: screwing up the init.amts field.")
  function fwd_sim(prev_sample_state::Array{Int64, 1}, prop_sample_params::Array{Float64, 1}, T_sim::Float64)
    for i in 1:length(unk_inds)
      wilk_cri.rxn_rates[unk_inds[i]] = prop_sample_params[i]
    end
    wilk_cri.init_amts = prev_sample_state
    return Chem_rxn_tools.gillespie(wilk_cri, T_sim, true)
  end
  MCS.fwd_sim = fwd_sim

#-----------------------------setup what and where to save-------------------------------
  ep.time_of_test = now()
  MCS.save_path = string(ep.save_folder, ep.time_of_test, "/")
  mkdir(MCS.save_path)

#-----------------------------Simulate the observations; plot; save plot-------------------------------
sim_results = Chem_rxn_tools.make_sim_data(t_obs, wilk_cri, obs_mol_name, noise_distribution)
  Chem_rxn_tools.plot_save_sim_data(MCS.save_path, sim_results, wilk_cri, mols_to_show)

metadata_path = joinpath(MCS.save_path, "everything_just_before_inference")
  save(metadata_path, exp_prefs, wilk_cri, sim_results)
  pMCMC_julia.MCS_save("before_inf_sample", MCS)

#-----------------------------do the inference-------------------------------
tic()
  pMCMC_julia.pMCMC!(sim_results.d_obs, sim_results.t_obs, MCS)
  ep.time_taken = toc()

#-----------------------------save/print the results------------------------------
  pMCMC_julia.MCS_save("samples_and_metadata", MCS)

  Prelim_experiments.print_exp_prefs(ep)
  println(string("By round, the acceptance rate was ", MCS.num_acc/(MCS.burnin_len + MCS.thin_len*num_samples_desired)))

end



#-----------------------------plot the posterior-------------------------------
ep = load(metadata_path)
string(MCS.save_path, "/samples_and_metadata")
  post_hists = pMCMC_julia.plot_save_marginals(MCS, MCS.save_path, unk_rates, sim_results.x_path[end])
  include("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/contour_bivariate_plot_maker.jl")
  post_contour = contour_plot_two_mols(MCS, MCS.save_path, unk_names[num_unks-1], unk_names[num_unks], wilk_cri, unk_rates, unk_names)
  include("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/bivariate_plot_maker.jl")
  post_biv = plot_save_two_mols(MCS, MCS.save_path, unk_names[num_unks-1], unk_names[num_unks], wilk_cri, unk_rates, unk_names)
  include("/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/project_specific/plot_biv_at_each_stage.jl")
  plot_biv_at_each_stage(MCS.save_path)

end

end
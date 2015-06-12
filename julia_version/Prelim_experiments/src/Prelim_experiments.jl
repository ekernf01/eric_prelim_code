module Prelim_experiments
  using Chem_rxn_tools
  using Distributions
  using pMCMC_julia
  using Dates
  using HDF5, JLD
  using Winston
  using KernelDensity


type Exp_prefs
  num_samples_desired::Int64
  bw_max::Float64
  bw_min::Float64
  bandwidth_multiplier::Float64

  true_noise_sd::Float64
  ass_noise_sd::Float64

  cri_source::String
  verbose::Bool

  obs_mol_name::String
  obs_mol_ind::Int64
  mols_to_show::Array{String,1}

  unk_names::Array{String,1}
  unk_rates::Array{Float64,1}
  unk_inds::Array{Int64,1}
  num_unks::Int64

  t_interval::Float64
  num_intervals::Int64
  save_folder::String

  time_of_test::DateTime
  time_taken::Float64
end

function Exp_prefs(num_samples_desired,
  bw_max,
  bw_min,
  bandwidth_multiplier,
  true_noise_sd,
  ass_noise_sd,
  cri_source,
  obs_mol_name,
  mols_to_show,
  unk_names,
  t_interval,
  num_intervals,
  save_folder)
  return (Exp_prefs(num_samples_desired,
          bw_max,
          bw_min,
          bandwidth_multiplier,

          true_noise_sd,
          ass_noise_sd,

          cri_source,
          false,       #verbose default

          obs_mol_name,
          0,          #obs_mol_ind to be filled in later
          mols_to_show,
          unk_names,
          Float64[], #unk_rates
          Int64[], #unk_inds
          length(unk_names), #num_unks

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
  println("verbose: ", ep.verbose)
  println("obs_mol_name: ", ep.obs_mol_name)
  println("unk_names: ", ep.unk_names)
  println("t_interval: ", ep.t_interval)
  println("num_intervals: ", ep.num_intervals)
  println("save_folder: ", ep.save_folder)
  println("time_of_test: ", ep.time_of_test)
  println("time_taken: ", ep.time_taken)
end

include("plot_biv_at_each_stage.jl")
include("plot_state_CI.jl")
include("start_off_test.jl")
include("pick_up_test.jl")

function make_all_plots(save_path)
  #-----------------------------plot the posterior-------------------------------
  wd = pwd()
  cd()
    cd(save_path)
    ep = load("everything_just_before_inference", "ep")
    MCS = pMCMC_julia.MCS_load("after_inf_sample")
    plot_biv_at_each_stage(save_path)
    #space to add time-series and histogram plots later
  cd()
  cd(wd)
  #plot_save_trajec_CIs(save_path)
end

end

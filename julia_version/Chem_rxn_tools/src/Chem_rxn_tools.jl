module Chem_rxn_tools
println("Apologies: Chem_rxn_tools takes a little while to load.")

#-------------------------------------Data structures-------------------------------------
#-------------------------------------CRI-------------------------------------
type Chem_rxn_info
  #Info about relevant chemicals
  species_labels::Array{String, 1}
  init_amts::Array{Int64, 1}
  num_species::Int64

  #Stoichiometry and reactions
  sto_mat_nonzero_inds::Array{Array{Int64, 1}}
  rxn_entry_mat_nonzero_inds::Array{Array{Int64, 1}}
  sto_mat::Matrix{Int64}
  rxn_entry_mat::Matrix{Int64}
  rxn_labels::Array{String,1}
  num_rxns::Int64
  rxn_rates::Array{Float64, 1}
  rxns_written_out::Array{String,1}

  #Data internal to the SMBL representation but otherwise extraneous for this module
  rxn_pos_in_SBML_file::Array{Int64,1}
  SBML_par_names::Array{String, 1}
  SBML_par_vals::Array{Float64, 1}
end

#constructor to make an empty chem_rxn_info
function Chem_rxn_info()
  println("Warning: always run chem_rxn_data_check!() after populating a Chem_rxn_info object.")
  return (Chem_rxn_info(
  String[],                            #species_labels
  Int64[],                             #init_amts
  0,                                   #num_species

  Array{Int64}[],                      #sto_mat_nonzero_inds
  Array{Int64}[],                      #rxn_entry_mat_nonzero_inds
  zeros(Int64, 0, 0),                  #sto_mat
  zeros(Int64, 0, 0),                  #rxn_entry_mat

  String[],                            #rxn_labels
  0,                                   #num_rxns
  Float64[],                           #rxn_rates
  String[],                            #rxns_written_out

  Int64[],                             #rxn_pos_in_SBML_file
  String[],                            #SBML_par_names
  Float64[]                            #SBML_par_vals
  ))
end

#constructor to mimic old constructor
function Chem_rxn_info(
  species_labels,
  init_amts,
  num_species,

  sto_mat,
  rxn_entry_mat,

  rxn_labels,
  num_rxns,
  rxn_rates,
  rxns_written_out,

  rxn_pos_in_SBML_file,
  SBML_par_names,
  SBML_par_vals)

  println("Warning: always run chem_rxn_data_check!() after populating a Chem_rxn_info object.")
  return (Chem_rxn_info(


  species_labels,
  init_amts,
  num_species,

  Array{Int64}[],                      #sto_mat_nonzero_inds
  Array{Int64}[],                      #rxn_entry_mat_nonzero_inds
  sto_mat,
  rxn_entry_mat,

  rxn_labels,
  num_rxns,
  rxn_rates,
  rxns_written_out,

  rxn_pos_in_SBML_file,
  SBML_par_names,
  SBML_par_vals))
end

#-------------------------------------Chem_sim_result-------------------------------------
type Chem_sim_result
  x_path::Array{Array{Int64, 1}, 1}
  current_x::Array{Int64, 1}
  num_rxns_occ::Int64
  rxn_types::Array{Int64, 1}
  rxn_times::Array{Float64, 1}
  t_spent::Float64
  x_obs::Array{Array{Int64, 1}, 1} #Noiseless verion of all values at observation times.
  d_obs::Array{Float64, 1} #Noisy verion of observed coordinate at observation times.
  t_obs::Array{Float64, 1} #observation times.
  obs_mol_name::String
  obs_mol_ind::Int64
end



using Distributions
using Winston

#-------------------------------------functions-------------------------------------
include("chem_rxn_data_check.jl")
include("get_chem_indices.jl")
include("get_rate_indices.jl")
include("get_par_info.jl")
include("get_species_info.jl")
include("get_rxn_info.jl")
include("get_stoich_info.jl")
include("get_rate_info.jl")
include("SBML_read.jl")
include("gillespie.jl")
include("gillespie_tester.jl")
include("make_cri_graphic.jl")
include("make_sim_data.jl")
include("plot_save_sim_data.jl")

#Produces a demo rxn system with specified number of molecules.
#slow rxn rates with decay 100x slower than immigration and init vals of 50.
#Interaction: remaining molecules catalyze creation of mol1.
function make_demo_cri_v1(num_species::Int64)

    demo_cri = Chem_rxn_info(
    [string("mol",i ) for i in 1:num_species],           #species_labels
    50*ones(Int64, num_species),                         #init_amts
    num_species,                                         #num_species

    hcat(eye(Int64,num_species, num_species),
     -eye(Int64,num_species, num_species),
     ones(Int64,num_species, 1)),                         #sto_mat
    hcat(zeros(Int64,num_species, num_species),
     eye(Int64,num_species, num_species),
     eye(Int64,num_species, 1)),                        #rxn_entry_mat

    vcat([string("prod_mol",i ) for i in 1:num_species],
    [string("decay_mol",i ) for i in 1:num_species],
         ["mol1_boost"]),                                #rxn_labels

    2*num_species,                                       #num_rxns
    vcat([0.01 for i in 1:num_species],
    [0.0001 for i in 1:num_species],
    [0.001]),                                            #rxn_rates
    vcat([string("\phi -> ", "mol",i) for i in 1:num_species],
    [string("mol",i, "-> \phi" ) for i in 1:num_species],
    ["mol1 splits into others"]),                        #rxns_written_out

    Int64[],                                             #rxn_pos_in_SBML_file
    String[],                                            #SBML_par_names
    Float64[]                                            #SBML_par_vals
    )
  chem_rxn_data_check!(demo_cri)
  return demo_cri
end

#Produces a demo rxn system with specified number of molecules.
#slow rxn rates with decay 100x slower than immigration and init vals of 50.
#Interaction: remaining molecules catalyze creation of mol1.
function make_demo_cri_v2()
    demo_cri = Chem_rxn_info(
      ["mol1"],                                              #species_labels
      500*ones(Int64, 1),                                   #init_amts
      1,                                                   #num_species

      reshape(Int64[-1, 1], 1, 2),                         #sto_mat
      eye(Int64, 1, 2),                                    #rxn_entry_mat
      ["decay", "prod"],                                   #rxn_labels
      2,                                                   #num_rxns
      [0.0001, 0.01],                                      #rxn_rates
      ["decay", "prod"],                                   #rxns_written_out

      Int64[],                                             #rxn_pos_in_SBML_file
      String[],                                            #SBML_par_names
      Float64[]                                            #SBML_par_vals
    )
  chem_rxn_data_check!(demo_cri)
  return demo_cri
end
end

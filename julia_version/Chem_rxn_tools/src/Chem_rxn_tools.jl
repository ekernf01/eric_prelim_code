module Chem_rxn_tools
# Module containing tools for:
# --storing info about chemical reactions
# --simulating chemical reactions via a markov jump process (implementation of the gillespie algorithm)
# --reading chemical reaction systems from SBML shorthand (SBML is "systems biology markup language")
# Example code:
# SBML_file = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/chem_rxn_tools/wilkinson_rxns_SBML_shorthand.txt"
# using Chem_rxn_tools
# wilk_cri = SBML_read(SBML_file)

  #-------------------------------------Data structures-------------------------------------
  type Chem_rxn_info
    #Info about relevant chemicals
    species_labels::Array{String, 1}
    init_amts::Array{Int64, 1}
    num_species::Int64

    #Stoichiometry and reactions
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
  Chem_rxn_info() = Chem_rxn_info(
    String[],                            #species_labels
    Int64[],                             #init_amts
    0,                                   #num_species

    zeros(Int64, 0, 0),                  #sto_mat
    zeros(Int64, 0, 0),                  #rxn_entry_mat
    String[],                            #rxn_labels
    0,                                   #num_rxns
    Float64[],                           #rxn_rates
    String[],                            #rxns_written_out

    [],                                  #rxn_pos_in_SBML_file
    String[],                            #SBML_par_names
    Int64[]                              #rxn_labels
    )

  immutable Chem_sim_result
    x_path::Array{Array{Int64, 1}, 1}
    current_x::Array{Int64, 1}
    num_rxns_occ::Int64
    rxn_types::Array{Int64, 1}
    rxn_times::Array{Float64, 1}
    t_spent::Float64
  end



  using Distributions
  using Winston

  #-------------------------------------functions-------------------------------------
  include("chem_rxn_data_check.jl")
  include("get_chem_indices.jl")
  include("get_par_info.jl")
  include("get_species_info.jl")
  include("get_rxn_info.jl")
  include("get_stoich_info.jl")
  include("get_rate_info.jl")
  include("SBML_read.jl")
  include("gillespie.jl")
  include("gillespie_tester.jl")
  include("make_cri_graphic.jl")
end


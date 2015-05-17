
using pMCMC_julia
pMCMC(d_obs, t_obs, current_MCMC_state::MCMC_state, emission_logden::Function, fwd_sim::Function)

using Distributions
using Winston
using ProfileView
using Dates
using HDF5, JLD
using Chem_rxn_tools
using pMCMC_julia


SBML_file = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/chem_rxn_tools/wilkinson_rxns_SBML_shorthand.txt"





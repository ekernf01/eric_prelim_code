using HDF5, JLD
using pMCMC_julia
MCS = pMCMC_julia.MCMC_state()

pMCMC_julia.MCS_save("hhhhgggg",MCS)
load_path = string(MCS.save_path, "hhhhgggg")
MCS = 0
MCS = pMCMC_julia.MCS_load(load_path)
MCS


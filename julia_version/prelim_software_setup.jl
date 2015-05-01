immutable Sample_state_and_params_type
  params::Array{Float64,2}
  state::Array{Float64,2}
end
#thetax = Sample_state_and_params_type([1 1;1 1], [1 1;1 1])
#thetax.params

using Distributions
include("gillespie.jl")
include("pMCMC.jl")


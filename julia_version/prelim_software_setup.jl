immutable Sample_state_and_params_type
  params::Array{Float64,2}
  state::Array{Int64,2}
end
Sample_state_and_params_type(len, par_dim, state_dim) =
  Sample_state_and_params_type(zeros(Float64,len, par_dim), zeros(Int64,len, state_dim))
#thetax = Sample_state_and_params_type([1 1;1 1], [1 1;1 1])
#thetax.params

using Distributions
using Winston
include("gillespie.jl")
include("pMCMC.jl")


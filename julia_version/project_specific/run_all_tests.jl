println("Testing GFP")
include("run_big_test_GFP_obs.jl")
println("Done testing GFP. Testing Hag.")
include("run_big_test_Hag_obs.jl")
println("Done testing GFP & Hag. Testing SigD.")
include("run_big_test_SigD_obs.jl")
println("Done testing GFP & Hag & SigD. Testing simple GFP.")
include("run_big_test_GFP_simple.jl")



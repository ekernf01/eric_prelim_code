workspace()
  push!(LOAD_PATH, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/")
  using pMCMC_julia
  using Winston
  function eric_qq(x, y,eric_title="")
    p = plot(sort(x), sort(y), title=eric_title)
    the_max = maximum([maximum(x),maximum(y)])
    the_min = minimum([minimum(x),minimum(y)])
    the_range = [the_min, the_max]
    add(p, Curve(the_range,the_range))
    return(p)
  end

samp_size = 10_000
  tic()
  (geweke_fwd_samp, geweke_mc_samp) = pMCMC_julia.simple_geweke_test(samp_size)
  time_taken = toc()
  println("Time taken: ", time_taken)

param_qq = eric_qq(geweke_fwd_samp[1][:], geweke_mc_samp[1][:], "Parameters") #params
state_qq = eric_qq(geweke_fwd_samp[2][:], geweke_mc_samp[2][:], "States") #states
obs_qq = eric_qq(geweke_fwd_samp[3][:], geweke_mc_samp[3][:], "Observations") #obs

savedir = "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/pMCMC_julia/tests"
  if(!isdir(savedir))
    mkdir(savedir)
  end
  savename = joinpath(savedir, "params.png")
  savefig(param_qq,savename)
  savename = joinpath(savedir, "state.png")
  savefig(state_qq,savename)
  savename = joinpath(savedir, "obs.png")
  savefig(obs_qq,savename)

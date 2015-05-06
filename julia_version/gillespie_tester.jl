include("gillespie.jl")

#Working example that calls gillespie
function gillespie_basic_test()
  init_x = Int [0,0]
  sto_mat = ones(Int,2,2)
  rxn_entry_mat = zeros(Int,2,2)
  rxn_rates = ones(2)
  T_sim = 5
  inside_sampler = false
  x_path, current_x, num_rxns_occ, rxn_types,rxn_times, t_spent = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim, inside_sampler)
  gillespie_input_ok(init_x, sto_mat, rxn_entry_mat, rxn_rates, 5.0)

  if rxn_times[1] != 0
    println("rxn_times did not have first element zero.")
  end
  if x_path[1] != init_x
    println("x_path did not have first element init_x.")
  end
  return 0
end

gillespie_basic_test()


# Unit tests should:
# --feed it garbage: strings, negative numbers and fractions where only positive numbers or integers should be, matrices where it needs vectors, stuff that's the wrong length
# --Run it many times on a simple reaction model and check to see it converge to the true stochastic mean
# --Run it and compare to Wilkinson's Figure 2
# --Make sure:
#      t_spent is always above T_sim
#      rxn_times increases but never hits T_sim
#      rxn_types add up properly
#      t_spent is always above t_sim
#      all output lines up with num_rxns_occ
# ——————————————————————————————————————————————

function gillespie_garbage_initx_test()
  init_x = 'a'
  sto_mat = ones(Int64,2,2)
  rxn_entry_mat = zeros(Int64,2,2)
  rxn_rates = ones(2)
  T_sim = 5
  inside_sampler = false
  all_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim, inside_sampler)
  return 0
end

function gillespie_garbage_stomat_test()
  init_x = Int64 [0,0]
  sto_mat = -ones(Int64,2,2)
  rxn_entry_mat = zeros(Int64,2,2)
  rxn_rates = ones(2)
  T_sim = 5
  inside_sampler = false
  all_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim, inside_sampler)
  return 0
end

function gillespie_garbage_rxn_entry_test()
  init_x = Int64 [0,0]
  sto_mat = ones(Int64,2,2)
  rxn_entry_mat = zeros(Int64,2,3)
  rxn_rates = ones(2)
  T_sim = 5
  inside_sampler = false
  all_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim, inside_sampler)
  print(all_out)
end

function gillespie_garbage_rxn_rates_test()
  init_x = Int64 [0,0]
  sto_mat = ones(Int64,2,2)
  rxn_entry_mat = zeros(Int64,2,2)
  rxn_rates = -ones(2)
  T_sim = 5
  inside_sampler = false
  all_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim, inside_sampler)
  print(all_out)
end

function gillespie_garbage_rxn_rates_2_test()
  init_x = Int64 [0,0]
  sto_mat = ones(Int64,2,2)
  rxn_entry_mat = zeros(Int64,2,2)
  rxn_rates = ones(3)
  T_sim = 5
  inside_sampler = false
  all_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim, inside_sampler)
  return 0
end

function gillespie_garbage_T_sim_test()
  init_x = Int64 [0,0]
  sto_mat = ones(Int64,2,2)
  rxn_entry_mat = zeros(Int64,2,2)
  rxn_rates = ones(3)
  T_sim = -5
  inside_sampler = false
  all_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim, inside_sampler)
  return 0
end

using Winston

function gillespie_stoch_mean_test(num_runs, k, n_0)
  #Simplest reaction: exponential decay
  T_sim = 50
  tgrid = [0:0.1:T_sim]
  init_x = [100]
  sto_mat = -ones(Int64,1,1)
  rxn_entry_mat = ones(Int64,1,1)
  rxn_rates = 0.01*ones(Float64, 1)
  inside_sampler = false

  #initialize plot
  decay_plot = FramedPlot(
           title="Eric Tests Gillespie Implementation",
           xlabel="Rxn Time",
           ylabel="Remaining Molecules")

  #plot many runs of the gillespie algo
  for which_run = [1:num_runs]
    x_path, current_x, num_rxns_occ, rxn_types,rxn_times, t_spent = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim, inside_sampler)
    add(decay_plot,Curve(rxn_times, [x_path...]))
  end

  #plot the analytically computable stochastic mean
  tgrid = [0:0.1:T_sim]
  add(decay_plot,Curve(tgrid, init_x.*exp(sto_mat[1,1].*rxn_rates[1,1].*tgrid), "color","red"))
  return decay_plot
end
decay_plot = gillespie_stoch_mean_test(20, 0.01*ones(Float64, 1), 100)
savefig(decay_plot, "/Users/EricKernfeld/Desktop/Spring_2015/518/eric_prelim_code/julia_version/gillespie_test_figure.pdf")



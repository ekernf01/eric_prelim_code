
include("gillespie.jl")

#Working example that calls gillespie
function gillespie_basic_test()
  init_x = Int64 [0,0]
  sto_mat = ones(Int64,2,2)
  rxn_entry_mat = zeros(Int64,2,2)
  rxn_rates = ones(2)
  T_sim = 5
  all_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim)
  return 0
end

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
  all_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim)
  return 0
end

function gillespie_garbage_stomat_test()
  init_x = Int64 [0,0]
  sto_mat = -ones(Int64,2,2)
  rxn_entry_mat = zeros(Int64,2,2)
  rxn_rates = ones(2)
  T_sim = 5
  all_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim)
  return 0
end

function gillespie_garbage_rxn_entry_test()
  init_x = Int64 [0,0]
  sto_mat = ones(Int64,2,2)
  rxn_entry_mat = zeros(Int64,2,3)
  rxn_rates = ones(2)
  T_sim = 5
  all_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim)
  print(all_out)
end

function gillespie_garbage_rxn_rates_test()
  init_x = Int64 [0,0]
  sto_mat = ones(Int64,2,2)
  rxn_entry_mat = zeros(Int64,2,2)
  rxn_rates = -ones(2)
  T_sim = 5
  all_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim)
  print(all_out)
end

function gillespie_garbage_rxn_rates_2_test()
  init_x = Int64 [0,0]
  sto_mat = ones(Int64,2,2)
  rxn_entry_mat = zeros(Int64,2,2)
  rxn_rates = ones(3)
  T_sim = 5
  all_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim)
  return 0
end

function gillespie_garbage_T_sim_test()
  init_x = Int64 [0,0]
  sto_mat = ones(Int64,2,2)
  rxn_entry_mat = zeros(Int64,2,2)
  rxn_rates = ones(3)
  T_sim = -5
  all_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim)
  return 0
end

using Gadfly
#function gillespie_stoch_mean_test(num_runs, k, n_0)
  #Simplest reaction: exponential decay
  T_sim = 50
  tgrid = [0:0.1:T_sim]
  #for which_run = [1:num_runs]
    init_x = 100
    sto_mat = -ones(Int64,1,1)
    rxn_entry_mat = ones(Int64,1,1)
    rxn_rates = 0.01*ones(Float64, 1)
    gillespie_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim)
  #end
  plot(x=gillespie_out[4], y=collect([1:length(gillespie_out[3])]))
  #return gillespie_out
#end
#gillespie_stoch_mean_test(1, 0.01*ones(Float64, 1), 100)
  plot(x=gillespie_out[4], y=collect([1:length(gillespie_out[3])]))


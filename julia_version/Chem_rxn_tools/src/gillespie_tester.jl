#Working example that calls gillespie
function gillespie_basic_test()
  cri = Chem_rxn_info()
  cri.init_amts = Int[0,0]
  cri.sto_mat = ones(Int,2,2)
  cri.rxn_entry_mat = zeros(Int,2,2)
  cri.rxn_rates = ones(2)
  chem_rxn_data_check!(cri)

  T_sim = 5.0
  inside_sampler = false
  demo_sim_results = gillespie(cri, T_sim, inside_sampler)

  if demo_sim_results.rxn_times[1] != 0
    println("rxn_times did not have first element zero.")
  end
  if demo_sim_results.x_path[1] != cri.init_amts
    println("x_path did not have first element init_x.")
  end
  return 0
end


using Winston
function gillespie_stoch_mean_test(num_runs, k, n_0)
  #Simplest reaction: exponential decay
  cri = Chem_rxn_info()

  cri.init_amts = [100]
  cri.sto_mat = -ones(Int64,1,1)
  cri.rxn_entry_mat = ones(Int64,1,1)
  cri.rxn_rates = 0.01*ones(Float64, 1)
  cri.num_rxns = 1
  cri.num_species = 1
  chem_rxn_data_check!(cri)

  T_sim = 50.0
  tgrid = [0:0.1:T_sim]
  inside_sampler = false

  #initialize plot
  decay_plot = FramedPlot(
           title="Eric Tests Gillespie Implementation",
           xlabel="Rxn Time",
           ylabel="Remaining Molecules")

  #plot many runs of the gillespie algo
  for which_run = [1:num_runs]
    sim_results_temp = gillespie(cri, T_sim, inside_sampler)
    sim_results_temp.num_rxns_occ
    add(decay_plot,Curve(sim_results_temp.rxn_times, [sim_results_temp.x_path...]))
  end

  #plot the analytically computable stochastic mean
  tgrid = [0:0.1:T_sim]
  add(decay_plot,Curve(tgrid, cri.init_amts.*exp(cri.sto_mat[1,1].*cri.rxn_rates[1,1].*tgrid), "color","red"))
  return decay_plot
end

function gillespie_stoch_mean_test_sampler_version(num_runs, k, n_0)
  #Simplest reaction: exponential decay
  cri = Chem_rxn_info()
  cri.init_amts = [100]
  cri.sto_mat = -ones(Int64,1,1)
  cri.rxn_entry_mat = ones(Int64,1,1)
  cri.rxn_rates = 0.01*ones(Float64, 1)
  cri.num_rxns = 1
  cri.num_species = 1
  T_sim = 50.0
  tgrid = [0:0.1:T_sim]
  inside_sampler = true
  chem_rxn_data_check!(cri)

  #initialize plot
  decay_plot = FramedPlot(
           title="Eric Tests Gillespie Implementation for the Version Run Inside the Sampler",
           xlabel="Rxn Time",
           ylabel="Remaining Molecules")

  #plot many runs of the gillespie algo
  mean_final = 0
  for which_run = [1:num_runs]
    current_x = gillespie(cri, T_sim, inside_sampler)
    add(decay_plot,Points([0,T_sim], [cri.init_amts,current_x]))
    mean_final = mean_final + current_x
  end
  mean_final = mean_final/num_runs

  #plot the analytically computable stochastic mean
  tgrid = [0,T_sim]
  add(decay_plot,Points(tgrid, cri.init_amts.*exp( cri.sto_mat[1,1].* cri.rxn_rates[1,1].*tgrid), "color","red"))
  add(decay_plot,Points(T_sim, mean_final, "color","blue"))
  return decay_plot
end



# gillespie_stoch_mean_test_sampler_version(100, 0.1, 50)
# gillespie_stoch_mean_test(100, 0.1, 50)
# gillespie_basic_test()


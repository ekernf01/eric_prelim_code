init_x = Int64 [0,0]
sto_mat = ones(Int64,2,2)
rxn_entry_mat = zeros(Int64,2,2)
rxn_rates = ones(2)
T_sim = 5
all_out = gillespie(init_x, sto_mat, rxn_entry_mat, rxn_rates, T_sim)
print(all_out)
# Unit tests should:
# --feed it garbage: strings, negative numbers and fractions where only positive numbers or integers should be, matrices where it needs vectors, stuff that's the wrong length
# --Run it many times on a simple reaction model and check to see it converge to the true stochastic mean
# --Make sure:
#      t_spent is always above T_sim
#      rxn_times increases but never hits T_sim
#      rxn_types add up properly
#      t_spent is always above t_sim
#      all output lines up with num_rxns_occ
# ——————————————————————————————————————————————

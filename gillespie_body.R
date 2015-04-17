# ——————————————————————————————————————————————
# Forward simulation 
# ——————————————————————————————————————————————
# Given:
# 
# init_x, the particle counts
# sto_mat, the stoichiometry matrix 
# rxn_rates, reaction rates, MEASURED IN INTENSITY PER SECOND
# T_sim, the time over which to simulate, MEASURED IN SECONDS
# 
# It verifies: 
# 
# init_x is nonnegative integers
# length of init_x matches num rows of sto_mat
# length of rxn_rates matches num cols of sto_mat
# T_sim is nonnegative
# 
# Then, it runs the Gillespie algorithm to produce:
# 
# rxn_time, vector of times at which reactions occur (for convenience, rxn_time[0] should be 0)
# rxn_type, list of reactions that occured (integers that index sto_mat)
# x_path, a list of vectors with zero element init_x and i element showing the molecule counts between rxn_time[i] and rxn_time[i].
# 
# ——————————————————————————————————————————————
gillespie_body <- cxxfunction(
                includes='#include <RcppArmadilloExtensions/sample.h>',
                signature(init_x_r="numeric", sto_mat_r="numeric", 
                          rxn_rates_r="numeric", T_sim_r="numeric"),
                plugin="Rcpp", 
                body= '

  NumericVector init_x_c(init_x_r), sto_mat_c(sto_mat_r);
  NumericVector rxn_rates_c(rxn_rates_r), T_sim_c(T_sim_r);
  
  NumericVector rxn_time;
  NumericVector rxn_type;
  NumericVector x_path;

  NumericVector alpha;
  double alpha_tot = 1;
  int num_rxn_types = rxn_rates_c.size(); 
  NumericVector rxn_type_temp_indic(num_rxn_types);
  int num_rxn_events = 0;
  int num_mols = init_x_c.size();
  double t_i = 0;
  double T_tot = T_sim_c[0];
  rxn_time[0] = 0;

  RNGScope scope;

  while(t_i<T_tot){
    for(int which_mol=0; which_mol<num_mols; which_mol++){
      alpha[which_mol] = 1;
    }
    t_i = t_i + rexp(1)[0];
    if(t_i<T_tot){
      rxn_type_temp_indic = RcppArmadillo::sample(num_rxn_types, 1, FALSE, alpha);
      for(int i; i<num_rxn_types; i++){
        if(rxn_type_temp_indic){rxn_type = i;}
      }
      num_rxn_events++;
      rxn_time[num_rxn_events] = t_i;
    }
    
  }

  return init_x_c;
')
print(gillespie_body(1,1,1,1))

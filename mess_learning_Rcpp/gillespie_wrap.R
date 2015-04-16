#wrapper for gillespie_guts_c function in file gillespie_guts_c.cpp
Rcpp::sourceCpp('~/Dropbox/518/code/gillespie_guts_c.cpp')
gillespie_wrap <- function(init_x, sto_mat,
                           rxn_rates,  T_sim) {
  #To do: error checking here
  out_df <- rcpp("gillespie_guts_c", 
                 init_x, sto_mat,
                 rxn_rates,  T_sim)
  return(out_df)
}

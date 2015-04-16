#include <Rcpp.h>
#include <R.h>
#include <Rmath.h>
using namespace Rcpp;

// Below is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar)

// For more on using Rcpp click the Help button on the editor toolbar
//
//——————————————————————————————————————————————
//Forward simulation 
//(name ideas: gillespie.c, fwd_sim.c)
//——————————————————————————————————————————————
//Given:
//
//init_x, the particle counts
//sto_mat, the stoichiometry matrix 
//c, reaction rates, MEASURED IN INTENSITY PER SECOND
//T_sim, the time over which to simulate, MEASURED IN SECONDS
//
//It verifies: 
//
//init_x is nonnegative integers
//length of init_x matches num rows of sto_mat
//length of c matches num cols of sto_mat
//T_sim is nonnegative
//
//Then, it runs the Gillespie algorithm to produce:
//
//rxn_time, vector of times at which reactions occur (for convenience, rxn_time[0] should be 0)
//rxn_type, list of reactions that occured (integers that index sto_mat)
//x_path, a list of vectors with zero element init_x and i element showing the molecule counts between rxn_time[i] and rxn_time[i].
//
//——————————————————————————————————————————————
//

// [[Rcpp::export]]
SEXP gillespie_guts_c(SEXP init_x, 
                    SEXP sto_mat, 
                    SEXP rxn_rate, 
                    SEXP T_sim){
    GetRNGstate();//must precede all calls to R's RNG functions
    
    //Convert input from R to C
    Rcpp::NumericVector init_x_c(init_x);
    Rcpp::NumericMatrix sto_mat_c(sto_mat);
    Rcpp::NumericVector rxn_rate_c(rxn_rate);
    Rcpp::NumericVector T_sim_c(T_sim);
    
    //Generate output
    Rcpp::NumericVector rxn_time;
    Rcpp::NumericVector rxn_type;
    Rcpp::NumericMatrix x_path;
    
    int num_rxn = 0;
    Rcpp::NumericVector sim_time;
    //while(sim_time < T_sim){
        //replace this with the actual gillespie algo
        sim_time = sim_time + rexp(1);
        num_rxn++;
    //}
    PutRNGstate();//must follow all calls to R's RNG functions

    //Set up the output
    Rcpp::DataFrame out_df = Rcpp::DataFrame::create(
      Rcpp::Named("rxn_time")=rxn_time,
      Rcpp::Named("rxn_type")=rxn_type,
      Rcpp::Named("x_path")=x_path);
    
    Rcpp::List::create(
      Rcpp::Named("out_df")=out_df,
      Rcpp::Named("num_rxn")=num_rxn);   
                    
    return 0;
}
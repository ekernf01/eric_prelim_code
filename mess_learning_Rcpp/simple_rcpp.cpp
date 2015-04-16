#include <Rcpp.h>
 
  // [[Rcpp::export]]
extern "C" SEXP simple_rcpp(SEXP x){
    Rcpp::NumericVector x_c = x;
    return x;
}
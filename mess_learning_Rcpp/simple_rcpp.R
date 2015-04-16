#a simple Rcpp example from Dirk Eddelbuettel
ex1rcppSugar <- cxxfunction(
                signature(a="numeric", b="numeric"),
                plugin="Rcpp", 
                body= '
  NumericVector x(a), y(b);
  return x;
')
print(ex1rcppSugar(1,1))

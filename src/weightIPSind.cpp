#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

double innerWind(arma::mat& X, int j, int l, int npar) {
  double v = 1.0;
  
  int jj = j;
  int ll = l;
  
  for (int s = 0;  s < npar; s++ ){
    v = v * (X(jj,s) <= X(ll,s));
  }
  
  return v;
}


// [[Rcpp::export]]
arma::mat weightIPSind(arma::mat X) {
  int nobj = X.n_rows, npar = X.n_cols;
  
  arma::mat wind = eye<arma::mat>(nobj, nobj);
 
  for (int l = 0; l < nobj; l++){
    for (int j = 0; j < nobj; j++){
      wind(j,l) =  innerWind(X,j,l,npar);
    }
  }
      
  return (wind * trans(wind)/nobj);
}
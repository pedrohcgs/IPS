#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;



// [[Rcpp::export]]
arma::mat weightIPSexp(arma::mat X, String X_trans) {
  int nobj = X.n_rows, npar = X.n_cols;
  double v = 0;
  arma::mat wexp = eye<arma::mat>(nobj, nobj);;
  arma::mat xprob(nobj, npar);
  arma::mat vcinv(npar, npar);

  
  if (X_trans == "normal") {
    
    for (int l = 0;  l < npar; l++ ){
      xprob.col(l) = normcdf(( X.col(l) - mean(X.col(l)) ) / stddev(X.col(l)));
      
    }

    for (int j = 0; j < (nobj-1); j++){
      for (int l = (j+1); l < nobj ; l++){
        v = as_scalar(dot(xprob.row(j) - xprob.row(l), xprob.row(j) - xprob.row(l))) ;
        wexp(j,l) = exp(-0.5 * v);
        wexp(l,j) = wexp(j,l);
      }
    }
    
    
  } else if (X_trans == "arctan") {
    
    for (int l = 0;  l < npar; l++ ){
      xprob.col(l) = atan(( X.col(l) - v * mean(X.col(l)) ) / stddev(X.col(l)));
    }
    

    for (int j = 0; j < (nobj-1); j++){
      for (int l = (j+1); l <nobj ; l++){
        wexp(j,l) = exp(-0.5* as_scalar(dot(xprob.row(j) - xprob.row(l), xprob.row(j) - xprob.row(l))));
        wexp(l,j) = wexp(j,l);
      }
    }
    
  } else {
    vcinv = pinv(cov(X));
    

    for (int j = 0; j < (nobj-1); j++){
      for (int l = (j+1); l <nobj ; l++){
        wexp(j,l) = as_scalar((X.row(j) - X.row(l)) * vcinv * trans((X.row(j) - X.row(l))));
        wexp(j,l) = exp(-0.5 * wexp(j,l));
        wexp(l,j) = wexp(j,l);
      }
    }
    
  }
  return wexp;
}

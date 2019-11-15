#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
double objIPS(arma::vec b, arma::vec d, arma::mat X, arma::mat w, double  treated_flag, arma::vec whs) {
  int nobj = X.n_rows;
  double obj = 0.0;

  
  arma::vec psfit = X * b;
  psfit = 1.0/(1.0 + exp(-psfit));
  
  
  arma::vec h1 = (whs % d/psfit) / mean(whs % d/psfit) - 1.0;
  arma::vec h0 = (whs % (1.0 - d)/(1.0 - psfit))/mean(whs % (1.0 - d)/(1.0 - psfit)) - 1.0;
  
  obj = (1.0- treated_flag) * as_scalar(arma::trans(h1)* (w) * (h1))/(nobj*nobj) + 
          as_scalar(arma::trans(h0)* (w) * (h0))/(nobj*nobj);
  
  return obj;
}

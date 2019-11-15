#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec gradIPS(arma::vec b, arma::vec d, arma::mat& X, arma::mat& w, double treated_flag, arma::vec whs) {
  int nobj = X.n_rows, npar = X.n_cols;
  
  
  arma::vec psfit = X * b;
  psfit = 1.0/(1.0 + exp(-psfit));
 
  arma::vec w1 = (whs % d/psfit) / mean(whs % d/psfit);
  arma::vec w0 = (whs % (1.0 - d)/(1.0 - psfit)) / mean(whs % (1.0 - d)/(1.0 - psfit));
  
  // calculate h1.dot and h0.dot
  arma::mat h1dot(nobj, npar);
  arma::mat h0dot(nobj, npar);
  
  for (int j = 0;  j<npar; j++){
    h1dot.col(j) = - w1 % (1.0 - psfit) % X.col(j);
    h1dot.col(j) =  h1dot.col(j) - w1 * mean(h1dot.col(j));

    h0dot.col(j) =  w0 % psfit % X.col(j);
    h0dot.col(j) =  h0dot.col(j) - w0 * mean(h0dot.col(j));
  }

  arma::mat Qdot = ( (1.0 - treated_flag) * (arma::repmat(w1,1,nobj) -1.0) % w * h1dot + (arma::repmat(w0,1,nobj) -1.0) % w * h0dot)/nobj;
  arma::vec Qd = arma::trans(mean(Qdot, 0));
  return Qd;
}

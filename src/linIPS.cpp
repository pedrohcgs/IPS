#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat linIPS(arma::vec bhat, arma::vec d, arma::vec pshat, arma::mat& X, arma::mat& w, double treated_flag, arma::vec whs) {
  int nobj = X.n_rows, npar = X.n_cols;
  
  arma::mat Cinv = zeros<arma::mat>(npar,npar);
  arma::mat Cmat = zeros<arma::mat>(npar,npar);
  arma::mat linrep = zeros<arma::mat>(nobj,npar);
  
  arma::vec w1   = (whs % d/pshat) / mean(whs % d/pshat);
  arma::vec w0   = (whs % (1.0 - d)/(1.0 - pshat) )/ mean(whs % (1.0 - d)/(1.0 - pshat)) ;
  arma::mat h1dot = zeros<arma::mat>(nobj, npar);
  arma::mat h0dot= zeros<arma::mat>(nobj, npar);
  
  // calculate h11.dot
  for (int j = 0;  j< npar; j++){
    h1dot.col(j) = - w1 % (1.0 - pshat) % X.col(j);
    h1dot.col(j) =  h1dot.col(j) - w1 * mean(h1dot.col(j));

    h0dot.col(j) =  w0 % pshat % X.col(j);
    h0dot.col(j) =  h0dot.col(j) - w0 * mean(h0dot.col(j));
  }

  Cmat = ( (1.0 - treated_flag ) * trans(h1dot) * w * h1dot + trans(h0dot) * w * h0dot ) / nobj;
  bool flag = arma::inv(Cinv, Cmat);
  
  if (!flag) {
    Cinv = pinv(Cmat);
  } else {
    Cinv = arma::inv(Cmat);
  }
  
  linrep = - ( (1.0 - treated_flag) * (arma::repmat(w1,1,nobj) -1.0) % w * h1dot + (arma::repmat(w0,1,nobj) -1.0) % w * h0dot) * Cinv;
  return linrep;
}

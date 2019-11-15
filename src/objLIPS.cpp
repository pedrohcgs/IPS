#include <RcppArmadillo.h>
using namespace Rcpp;
using namespace arma;

// [[Rcpp::depends(RcppArmadillo)]]
// [[Rcpp::export]]
double objLIPS (arma::vec b, arma::vec d, arma::vec z, arma::mat& X, arma::mat& w,
                arma::vec whs) {
  int nobj = X.n_rows;
  double obj =  0;
  
  arma::vec ipsfit = X * b;
  ipsfit = 1.0/(1.0 + exp(-ipsfit));
  
  
  arma::vec hlte1(nobj, fill::zeros);
  arma::vec hlte0(nobj, fill::zeros);
  
  hlte1 =  (whs % d % ( z/ipsfit -  (1.0-z)/(1.0-ipsfit)) )/mean(whs % d % ( z/ipsfit -  (1.0-z)/(1.0-ipsfit)) )  -  ( whs % (  1.0 - (1.0 - d)%z/ipsfit - d%(1.0 - z)/(1.0 - ipsfit) ) ) / mean( whs % (  1.0 - (1.0 - d)%z/ipsfit - d%(1.0 - z)/(1.0 - ipsfit) ) );
  hlte0 =  (whs % (1.0-d) % ( (1.0-z)/(1.0-ipsfit) - z/ipsfit) )/mean(whs % (1.0-d) % ( (1.0-z)/(1.0-ipsfit) - z/ipsfit) )  -  ( whs % (  1.0 - (1.0 - d)%z/ipsfit - d%(1.0 - z)/(1.0 - ipsfit) ) )/ mean( whs % (  1.0 - (1.0 - d)%z/ipsfit - d%(1.0 - z)/(1.0 - ipsfit) ) );
  
  
  obj = as_scalar(trans(hlte1)* w * hlte1 + trans(hlte0)* w * hlte0)/(nobj*nobj); 
  
  
  return obj;
}
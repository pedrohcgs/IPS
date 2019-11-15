#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// [[Rcpp::export]]
arma::mat linLIPS(arma::vec bhat, arma::vec d, arma::vec z, arma::vec ipshat,
                  arma::mat& X, arma::mat& w, arma::vec whs) {
  int nobj = X.n_rows, npar = X.n_cols;
  
  arma::mat Cinv = zeros<arma::mat>(npar,npar);
  arma::mat Cmat = zeros<arma::mat>(npar,npar);
  arma::mat linrep = zeros<arma::mat>(nobj,npar);
  
  
  arma::vec hlte1(nobj, fill::zeros);
  arma::vec hlte0(nobj, fill::zeros);
  
  hlte1 =  (whs % d % ( z/ipshat -  (1.0-z)/(1.0-ipshat)) )/mean(whs % d % ( z/ipshat -  (1.0-z)/(1.0-ipshat)) )  -  ( whs % (  1.0 - (1.0 - d)%z/ipshat - d%(1.0 - z)/(1.0 - ipshat) ) ) / mean( whs % (  1.0 - (1.0 - d)%z/ipshat - d%(1.0 - z)/(1.0 - ipshat) ) );
  hlte0 =  (whs % (1.0-d) % ( (1.0-z)/(1.0-ipshat) - z/ipshat) )/mean(whs % (1.0-d) % ( (1.0-z)/(1.0-ipshat) - z/ipshat) )  -  ( whs % (  1.0 - (1.0 - d)%z/ipshat - d%(1.0 - z)/(1.0 - ipshat) ) )/ mean( whs % (  1.0 - (1.0 - d)%z/ipshat - d%(1.0 - z)/(1.0 - ipshat) ) );
  
  
  // calculate hlte dot
  arma::mat hltedot1(nobj, npar, fill::zeros);
  arma::mat hltedot0(nobj, npar, fill::zeros);
  arma::mat hltedotw1(nobj, npar, fill::zeros);
  arma::mat hltedotw(nobj, npar, fill::zeros);
  arma::mat hltedotw0(nobj, npar, fill::zeros);
  
  
  
  for (int j = 0;  j<npar; j++){
    hltedotw1.col(j) = - whs% d % (z % (1.0 - ipshat)/ipshat + (1.0 - z) % ipshat /(1.0 - ipshat)) /mean(whs % d % ( z/ipshat -  (1.0-z)/(1.0-ipshat)) ) %  X.col(j);
    hltedotw1.col(j) =  hltedotw1.col(j) +  whs % d % ( z/ipshat -  (1.0-z)/(1.0-ipshat))/mean(whs % d % ( z/ipshat -  (1.0-z)/(1.0-ipshat)) ) * mean(whs% d % (z % (1.0 - ipshat)/ipshat + (1.0 - z) % ipshat /(1.0 - ipshat)) /mean(whs % d % ( z/ipshat -  (1.0-z)/(1.0-ipshat)) ) %  X.col(j));
    
    hltedotw0.col(j) =  whs% (1.0-d) % (z % (1.0 - ipshat)/ipshat + (1.0 - z) % ipshat /(1.0 - ipshat)) /mean(whs % (1.0-d) % ( (1.0-z)/(1.0-ipshat) - z/ipshat) ) %  X.col(j);
    hltedotw0.col(j) =  hltedotw0.col(j) - whs % (1.0-d) % ( (1.0-z)/(1.0-ipshat) - z/ipshat)/mean(whs % (1.0-d) % ( (1.0-z)/(1.0-ipshat) - z/ipshat)) * mean( whs% (1.0-d) % (z % (1.0 - ipshat)/ipshat + (1.0 - z) % ipshat /(1.0 - ipshat)) /mean(whs % (1.0-d) % ( (1.0-z)/(1.0-ipshat) - z/ipshat) ) %  X.col(j));
    
    
    hltedotw.col(j) =   whs % ( (1.0 - d)% z % (1.0 - ipshat)/ipshat - d%(1.0 - z) % ipshat/(1.0 - ipshat) ) / mean(whs % (1.0 - (1.0 - d)%z/ipshat - d%(1.0 - z)/(1.0 - ipshat))) %  X.col(j); 
    hltedotw.col(j) =   hltedotw.col(j) -  whs%(1.0 - (1.0 - d)%z/ipshat - d%(1.0 - z)/(1.0 - ipshat))/mean(whs % (1.0 - (1.0 - d)%z/ipshat - d%(1.0 - z)/(1.0 - ipshat))) * mean(whs % ( (1.0 - d)% z % (1.0 - ipshat)/ipshat - d%(1.0 - z) % ipshat/(1.0 - ipshat) ) / mean(whs % (1.0 - (1.0 - d)%z/ipshat - d%(1.0 - z)/(1.0 - ipshat)))%  X.col(j));
    
    hltedot1.col(j) =  hltedotw1.col(j) - hltedotw.col(j);
    hltedot0.col(j) =  hltedotw0.col(j) - hltedotw.col(j);
    
  }
  
  
  
  
  Cmat = ( trans(hltedot1) * w *hltedot1 + trans(hltedot0) * w *hltedot0 )  / nobj;
  bool flag = arma::inv(Cinv, Cmat);
  
  if (!flag) {
    Cinv = pinv(Cmat);
    ///Rcout << "The Variance-Covariance Matrix is close to singular - proceed with caution!\n";
  } else {
    Cinv = arma::inv(Cmat);
  }
  
  linrep = -  ( arma::repmat(hlte1, 1, nobj)%w*hltedot1 + arma::repmat(hlte0, 1, nobj)%w*hltedot0) * Cinv;
  
  
  
  
  
  
  
  return linrep;
}
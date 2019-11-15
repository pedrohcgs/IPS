#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec gradLIPS(arma::vec b, arma::vec d, arma::vec z, arma::mat& X, arma::mat& w, 
                   arma::vec whs) {
  int nobj = X.n_rows, npar = X.n_cols;
  
  
  arma::vec ipsfit = X * b;
  ipsfit = 1.0/(1.0 + exp(-ipsfit));
  
  arma::vec hlte1(nobj, fill::zeros);
  arma::vec hlte0(nobj, fill::zeros);
  
  hlte1 =  (whs % d % ( z/ipsfit -  (1.0-z)/(1.0-ipsfit)) )/mean(whs % d % ( z/ipsfit -  (1.0-z)/(1.0-ipsfit)) )  -  ( whs % (  1.0 - (1.0 - d)%z/ipsfit - d%(1.0 - z)/(1.0 - ipsfit) ) ) / mean( whs % (  1.0 - (1.0 - d)%z/ipsfit - d%(1.0 - z)/(1.0 - ipsfit) ) );
  hlte0 =  (whs % (1.0-d) % ( (1.0-z)/(1.0-ipsfit) - z/ipsfit) )/mean(whs % (1.0-d) % ( (1.0-z)/(1.0-ipsfit) - z/ipsfit) )  -  ( whs % (  1.0 - (1.0 - d)%z/ipsfit - d%(1.0 - z)/(1.0 - ipsfit) ) )/ mean( whs % (  1.0 - (1.0 - d)%z/ipsfit - d%(1.0 - z)/(1.0 - ipsfit) ) );
  
  // calculate hlte dot
  arma::mat hltedot1(nobj, npar, fill::zeros);
  arma::mat hltedot0(nobj, npar, fill::zeros);
  arma::mat hltedotw1(nobj, npar, fill::zeros);
  arma::mat hltedotw(nobj, npar, fill::zeros);
  arma::mat hltedotw0(nobj, npar, fill::zeros);
  
  arma::mat Qdot(nobj, npar);
  arma::vec Qd(npar);
  
  
  for (int j = 0;  j<npar; j++){
    hltedotw1.col(j) = - whs% d % (z % (1.0 - ipsfit)/ipsfit + (1.0 - z) % ipsfit /(1.0 - ipsfit)) /mean(whs % d % ( z/ipsfit -  (1.0-z)/(1.0-ipsfit)) ) %  X.col(j);
    hltedotw1.col(j) =  hltedotw1.col(j) +  whs % d % ( z/ipsfit -  (1.0-z)/(1.0-ipsfit))/mean(whs % d % ( z/ipsfit -  (1.0-z)/(1.0-ipsfit)) ) * mean(whs% d % (z % (1.0 - ipsfit)/ipsfit + (1.0 - z) % ipsfit /(1.0 - ipsfit)) /mean(whs % d % ( z/ipsfit -  (1.0-z)/(1.0-ipsfit)) ) %  X.col(j));
    
    hltedotw0.col(j) =  whs% (1.0-d) % (z % (1.0 - ipsfit)/ipsfit + (1.0 - z) % ipsfit /(1.0 - ipsfit)) /mean(whs % (1.0-d) % ( (1.0-z)/(1.0-ipsfit) - z/ipsfit) ) %  X.col(j);
    hltedotw0.col(j) =  hltedotw0.col(j) - whs % (1.0-d) % ( (1.0-z)/(1.0-ipsfit) - z/ipsfit)/mean(whs % (1.0-d) % ( (1.0-z)/(1.0-ipsfit) - z/ipsfit)) * mean( whs% (1.0-d) % (z % (1.0 - ipsfit)/ipsfit + (1.0 - z) % ipsfit /(1.0 - ipsfit)) /mean(whs % (1.0-d) % ( (1.0-z)/(1.0-ipsfit) - z/ipsfit) ) %  X.col(j));
    
    
    hltedotw.col(j) =   whs % ( (1.0 - d)% z % (1.0 - ipsfit)/ipsfit - d%(1.0 - z) % ipsfit/(1.0 - ipsfit) ) / mean(whs % (1.0 - (1.0 - d)%z/ipsfit - d%(1.0 - z)/(1.0 - ipsfit))) %  X.col(j); 
    hltedotw.col(j) =   hltedotw.col(j) -  whs%(1.0 - (1.0 - d)%z/ipsfit - d%(1.0 - z)/(1.0 - ipsfit))/mean(whs % (1.0 - (1.0 - d)%z/ipsfit - d%(1.0 - z)/(1.0 - ipsfit))) * mean(whs % ( (1.0 - d)% z % (1.0 - ipsfit)/ipsfit - d%(1.0 - z) % ipsfit/(1.0 - ipsfit) ) / mean(whs % (1.0 - (1.0 - d)%z/ipsfit - d%(1.0 - z)/(1.0 - ipsfit)))%  X.col(j));
    
    hltedot1.col(j) =  hltedotw1.col(j) - hltedotw.col(j);
    hltedot0.col(j) =  hltedotw0.col(j) - hltedotw.col(j);
    
  }
  
  Qdot = (arma::repmat(hlte1,1,nobj)%w*hltedot1 + arma::repmat(hlte0,1,nobj)%w*hltedot0 ) /nobj;
  Qd = 2.0 * trans(mean(Qdot, 0));
  
  
  return Qd;
}
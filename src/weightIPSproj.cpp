#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;


// declare innerWproj
double innerWproj(arma::mat& X, int j, int l, int nobj) {
  double v = 0;
  double xjl = 0;
  double xjr = 0;
  double xlr = 0;
  double idxjl = 0;
  double arg = 0;
  double argbdd = 0;
  int jj = j;
  int ll = l;
  
  for (int r = 0;  r < nobj; r++ ){
    
    xjl = dot( (X.row(jj) - X.row(r)), (X.row(ll) - X.row(r)) );
    xjr = dot( (X.row(jj) - X.row(r)), (X.row(jj) - X.row(r)) );
    xlr = dot( (X.row(ll) - X.row(r)), (X.row(ll) - X.row(r)) );
    idxjl = dot( (X.row(jj) - X.row(ll)), (X.row(jj) - X.row(ll)) );
    
    if ((xlr == 0.0) && (xjr == 0.0)){
      v +=  2*M_PI;
    }
    else if (((xlr != 0.0) && (xjr != 0.0)) && (idxjl != 0.0)){
      arg = xjl/(sqrt(xjr)*sqrt(xlr));
      argbdd = std::min(std::max(arg, -1.0), 1.0);
      v += fabs(M_PI - std::acos(argbdd));
      arg = 0.0;
      argbdd = 0.0;
    }
    else {
      v +=  M_PI;
    }
    
  }
  
  return v;
}


// [[Rcpp::export]]
arma::mat weightIPSproj(arma::mat& X) {
  int nobj = X.n_rows;
  arma::mat wproj = zeros<arma::mat>(nobj, nobj);
 
  for (int j = 0;  j < nobj; j++ ){
    for (int l = (j); l <nobj ; l++){
      wproj(j,l) = innerWproj(X,j,l,nobj)/nobj;
      wproj(l,j) = wproj(j,l);
    }
  }
  
  return wproj;
}

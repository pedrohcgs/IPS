#include <Rcpp.h>
using namespace Rcpp;
#include <math.h>


// [[Rcpp::export]]
NumericMatrix weightIPSproj_uniq(const NumericMatrix& X, const NumericVector& wgt){
  int n_uniq = X.nrow();
  int k = X.ncol();
  NumericVector cum_wgt = cumsum(wgt);
  cum_wgt.push_front(0);
  int nobj = cum_wgt(n_uniq);
  
  double xjl;
  double xjr;
  double xlr;
  double idxjl;
  double arg;
  double argbdd;
  NumericMatrix omega(nobj, nobj);
  double out;
  
  xjl = 0.0;
  xjr = 0.0;
  xlr = 0.0;
  out = 0.0;
  idxjl = 0.0;
  arg = 0.0;
  argbdd = 0.0;
  
  for (int j = 0;  j<n_uniq; j++ ){
    for (int l = 0;  l<=j; l++ ){
      for (int r = 0;  r< n_uniq; r++ ){
        for (int s = 0; s < k; s++){
          xjl += (X(j,s) - X(r,s))*(X(l,s) - X(r,s));
          xjr += (X(j,s) - X(r,s))*(X(j,s) - X(r,s));
          xlr += (X(l,s) - X(r,s))*(X(l,s) - X(r,s));
          idxjl +=  (X(j,s) - X(l,s))*(X(j,s) - X(l,s));
        }
        
        if ((xlr == 0.0) && (xjr == 0.0)){
          out +=  2.0L * M_PI * wgt(r); }
        else if (((xlr != 0.0) && (xjr != 0.0)) && (idxjl != 0.0)){
          arg = (xjl / (sqrt(xjr)*sqrt(xlr)));
          if (arg > double(-1.0)){
            if (arg < double(1.0)){
              argbdd = arg;
            } else {
              argbdd = 1.0;
            }
          } else {
            argbdd = -1.0;
          }
          out += fabs(M_PI - acos(argbdd)) * wgt(r);
          arg = 0.0;
          argbdd = 0.0;
        }
        else {
          out +=  M_PI * wgt(r); }
        xjl = 0.0;
        xjr = 0.0;
        xlr = 0.0;
        idxjl = 0.0;
      }
      
      
      
      for (int bb = cum_wgt(j);  bb<= (cum_wgt(j+1)-1); bb++ ){
        for (int cc = cum_wgt(l);  cc<= (cum_wgt(l+1)-1); cc++ ){
        omega(bb, cc) = out;
        omega(cc, bb) = out;
      }
    }
      
      out = 0.0;
    }
    
  }
  
  
  return omega/nobj;
}












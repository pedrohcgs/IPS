#include <Rcpp.h>
using namespace Rcpp;
#include <math.h>


// [[Rcpp::export]]
NumericMatrix weightIPSproj(const NumericMatrix& X){
  int n = X.nrow();
  int k = X.ncol();
  double xjl;
  double xjr;
  double xlr;
  double idxjl;
  double arg;
  double argbdd;
  NumericMatrix omega(n,n);
  double out;
  
  xjl = 0.0;
  xjr = 0.0;
  xlr = 0.0;
  out = 0.0;
  idxjl = 0.0;
  arg = 0.0;
  argbdd = 0.0;
  
  for (int j = 0;  j<n; j++ ){
    for (int l = 0;  l<=j; l++ ){
      for (int r = 0;  r< n; r++ ){
        for (int s = 0; s < k; s++){
          xjl += (X(j,s) - X(r,s))*(X(l,s) - X(r,s));
          xjr += (X(j,s) - X(r,s))*(X(j,s) - X(r,s));
          xlr += (X(l,s) - X(r,s))*(X(l,s) - X(r,s));
          idxjl +=  (X(j,s) - X(l,s))*(X(j,s) - X(l,s));
        }
        
        if ((xlr == 0.0) && (xjr == 0.0)){
          out +=  2.0L * M_PI; }
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
          out += fabs(M_PI - acos(argbdd));
          arg = 0.0;
          argbdd = 0.0;
        }
        else {
          out +=  M_PI; }
        xjl = 0.0;
        xjr = 0.0;
        xlr = 0.0;
        idxjl = 0.0;
      }
      omega(j,l) = out;
      omega(l,j) = out;
      out = 0.0;
    }

  }
  
  
  return omega/n;
}












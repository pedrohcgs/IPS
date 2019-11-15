
// Code based on Adot, from the goffda R package

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
using namespace Rcpp;
using namespace arma;

// [[Rcpp::export]]
arma::vec weightIPSproj_vec(arma::mat X) {
  arma::uword n = X.n_rows;
  arma::mat inprod_full = arma::trimatu(X * X.t(), 0);
  arma::uvec ind_tri = arma::find(arma::trimatu(arma::ones(n, n), 0));
  arma::vec inprod = inprod_full.elem(ind_tri);
  arma::vec w_vec = arma::zeros(n * (n - 1) / 2 + 1);
  w_vec[0] = M_PI * (n + 1);

  for (arma::uword i = 2; i <= n; i++) {
    for (arma::uword j = 1; j < i; j++) {
      double sum_r = 0;
      for (arma::uword r = 1; r <= n; r++) {
        if ((i == r) | (j == r)) {
          sum_r += M_PI;
        } else {
          arma::uword aux_i = i * (i - 1) / 2;
          arma::uword aux_j = j * (j - 1) / 2;
          arma::uword aux_r = r * (r - 1) / 2;

          arma::uword ij = aux_i + j;
          arma::uword ii = aux_i + i;
          arma::uword jj = aux_j + j;
          arma::uword rr = aux_r + r;

          arma::uword ir = 0;
          if (i > r) {
            ir = aux_i + r;
            } else {

            ir = aux_r + i;

          }
        arma::uword rj = 0;
          if (r > j) {
            rj = aux_r + j;
          }
          else {
            rj = aux_j + r;
          }


          arma::uword jr = rj;
          double quo = inprod[ij - 1] - inprod[ir - 1] - inprod[rj - 1] + inprod[rr - 1];
          quo /=  std::sqrt((inprod[ii - 1] - 2 * inprod[ir - 1] + inprod[rr - 1]) *
            (inprod[jj - 1] - 2 * inprod[jr - 1] + inprod[rr - 1]));
          quo = std::min(std::max(quo, -1.0), 1.0);
          sum_r += std::acos(-quo);
        }
      }
      arma::uword ind_ij = (i - 1) * (i - 2) / 2 + j;
      w_vec[ind_ij] = sum_r;
    }
  }

  if (!w_vec.is_finite()) {
    warning("Non-finite values in weightIPSproj_vec");
  }
  return w_vec;
}

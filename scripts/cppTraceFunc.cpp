#include <RcppArmadillo.h>
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]

void showDims(const arma::mat& X){
  cerr << X.n_rows << " " << X.n_cols << endl;
}

// [[Rcpp::export]]
Rcpp::List bpFixedPointTrace(const arma::mat& bulk, const arma::mat& ref, int n_iter = 20, int window_size = 10) {

  int ncts = ref.n_cols;
  arma::mat ref_norm = ref.each_row() / arma::sum(ref, 0);
  double offset = arma::min(ref_norm.elem(find(ref_norm > 0)));
  ref_norm += offset;
  arma::mat refPnorm = ref_norm.each_col() / arma::sum(ref_norm, 1);

  arma::vec ppguess = arma::vec(ncts, arma::fill::ones) / ncts;
  arma::vec ppguess_pre = arma::vec(ncts);

  // Calculate the size of the log_likelihood vector based on window_size
  int ll_size = (n_iter % window_size == 0) ? n_iter / window_size : n_iter / window_size + 1;
  arma::vec log_likelihood = arma::vec(ll_size);
  int ll_index = 0;

  for (int i = 0; i < n_iter; ++i) {
    ppguess_pre = ppguess;
    arma::mat thisX = refPnorm.each_row() % ppguess.t();
    thisX = thisX.each_col() / arma::sum(thisX, 1);
    arma::mat thisXtot = thisX.each_col() % bulk;
    ppguess = arma::sum(thisXtot, 0).t();
    ppguess /= arma::sum(ppguess);

    // Use window_size to determine when to calculate log_likelihood
    if (i % window_size == 0 || i == n_iter - 1) {
      arma::mat log_p = ref.each_row() % ppguess.t();
      log_p = arma::pow(log_p, thisXtot);
      double offset2 = arma::min(log_p.elem(find(log_p > 0)));
      log_p += offset2;
      log_p = arma::log(log_p);
      log_likelihood[ll_index++] = arma::accu(log_p);
    }
  }

  return Rcpp::List::create(Rcpp::Named("pp") = ppguess, Rcpp::Named("log_likelihood") = log_likelihood);
}

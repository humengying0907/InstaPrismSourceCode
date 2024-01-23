#include <RcppArmadillo.h>
using namespace std;
// [[Rcpp::depends(RcppArmadillo)]]

//' @title bpFixedPoint function using Rcpp and Armadillo
//' @param bulk Numeric matrix
//' @param ref Numeric matrix
//' @param n_iter Number of iterations (default is 20)
//' @param ppguess Optional initial guess for proportion values
//' @return List with pp and perCell values
//' @export






void showDims(const arma::mat& X){
  cerr<<X.n_rows << "  " <<X.n_cols<<endl;
}

// [[Rcpp::export]]
Rcpp::List bpFixedPointTrace(const arma::mat& bulk, const arma::mat& ref, int n_iter = 20) {

  // Get number of columns in ref matrix
  int ncts = ref.n_cols;

  // Normalize column to sum to 1
  arma::mat ref_norm = ref.each_row() / arma::sum(ref, 0);

  // Add offset to make all elements positive
  double offset = arma::min(ref_norm.elem(find(ref_norm > 0)));
  ref_norm = ref_norm + offset;

  // Normalize rows to sum to 1
  arma::mat refPnorm = ref_norm.each_col()/ arma::sum(ref_norm, 1);

  //initial guess
  arma::vec ppguess = arma::vec(ncts, arma::fill::ones) / ncts;

  arma::vec ppguess_pre = arma::vec(ncts);
  arma::vec log_likelihood = arma::vec(n_iter);

  arma::mat thisX;
  arma::mat thisXtot;
  arma::mat log_p;

  // Iterate n_iter times
  for (int i = 0; i < n_iter; ++i) {

    ppguess_pre = ppguess;

    //multiply each row by current proportions
    thisX = refPnorm.each_row() % ppguess.t();

    // normalize each row to sum to 1
    thisX = thisX.each_col()/ arma::sum(thisX, 1);

    thisXtot = thisX.each_col() % bulk;

    ppguess = arma::sum(thisXtot, 0).t();

    ppguess = ppguess / arma::sum(ppguess);

    // calculate log likelihood
    log_p = ref.each_row() % ppguess.t();
    log_p = arma::pow(log_p, thisXtot);
    double offset2 = arma::min(log_p.elem(find(log_p > 0)));
    log_p = log_p + offset2;
    log_p = arma::log(log_p);
    log_likelihood[i] = arma::accu(log_p);

  }

  // Return result
  return Rcpp::List::create(Rcpp::Named("pp") = ppguess,
                            Rcpp::Named("log_likelihood") = log_likelihood);
}

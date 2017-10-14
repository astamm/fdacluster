#include "RcppArmadillo.h"
//[[Rcpp::depends(RcppArmadillo)]]


/// Main function
/**
 * This i s tha main function that will be exported to R by Rcpp.
 * For input and output look at R documentation.
 */

// [[Rcpp::export]]
Rcpp::List kmap(const arma::mat& x,
                const arma::cube& y,
                const arma::rowvec& seeds,
                const arma::uword n_clust,
                const std::string warping_method,
                const std::string center_method,
                const std::string similarity_method,
                const std::string optim_method,
                const arma::rowvec& warping_opt,
                const arma::rowvec& center_opt,
                const arma::rowvec& out_opt,
                const bool fence,
                const bool check_total_similarity,
                const bool show_iter,
                const bool comp_original_center,
                const arma::urowvec par_opt);

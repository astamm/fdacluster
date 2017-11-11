// Copyright (C) 2017 Alessandro Zito (zito.ales@gmail.com)
//
// This file is part of Fdakmapp.
//
// Fdakmapp is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Fdakmapp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Fdakmapp.  If not, see <http://www.gnu.org/licenses/>.

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

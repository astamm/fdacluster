// Copyright (C) 2017 Alessandro Zito (zito.ales@gmail.com)
//
// This file is part of fdacluster.
//
// fdacluster is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// fdacluster is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with fdacluster.  If not, see <http://www.gnu.org/licenses/>.

#include <RcppArmadillo.h>

/// Main function
/**
 * This is the main function that will be exported to R by Rcpp.
 * For input and output look at R documentation.
 */

// [[Rcpp::export]]
Rcpp::List kmap(
    const arma::mat &x,
    const arma::cube &y,
    const arma::urowvec &seeds,
    const arma::rowvec &warping_options,
    const unsigned int &n_clust,
    const unsigned int &maximum_number_of_iterations,
    const unsigned int &number_of_threads,
    const unsigned int &parallel_method,
    const double &distance_relative_tolerance,
    const double &center_args,
    const bool &cluster_on_phase,
    const bool &use_fence,
    const bool &check_total_dissimilarity,
    const bool &use_verbose,
    const bool &compute_overall_center,
    const std::string &warping_method,
    const std::string &center_method,
    const std::string &dissimilarity_method,
    const std::string &optimizer_method
);

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

#ifndef NEW_CENTERS_HPP
#define NEW_CENTERS_HPP

#include <RcppArmadillo.h>

#include "dissimilarity.h"
#include "center_methods.h"

/// Compute new centers
/**
 * It's a function that handle the computation of new centers.
 */

void newCenters(const arma::mat& x_reg,
                const arma::cube& y,
                const arma::rowvec& x_out,
                std::shared_ptr<Dissimilarity>& dissim,
                std::shared_ptr<CenterMethod>& cen,
                const arma::urowvec& par_opt,
                cube& templates,
                const urowvec& ict,
                const urowvec& labels,
                const bool show_iter
               );

#endif

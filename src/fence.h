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

#ifndef FENCE_HPP
#define FENCE_HPP

#include <RcppArmadillo.h>
#include "warping.h"
#include "optimizer.h"

/// Remove warping outliers
/** It's an optional check that can be activated by input. After each computation
 *  of best warping if the computed parameters are outliers they are recomputed with
 *  stricter bounds. It\' s computational less expensive way to have lighter warping
 *  insted of decrease input bound (for example: max_shift and max_dilation).
 */
void iterativeFence(
    arma::mat parameters,
    const arma::uword iter,
    arma::urowvec& labels,
    arma::rowvec& index,
    std::shared_ptr<WarpingFunction>& warping,
    std::shared_ptr<OptimizerMethod>& optimizer,
    const arma::cube& templates,
    const arma::mat& x_reg,
    const arma::cube& y,
    const arma::rowvec& x_out,
    std::shared_ptr<Dissimilarity>& dissim,
    const arma::urowvec& ict,
    const bool show_iter);
#endif

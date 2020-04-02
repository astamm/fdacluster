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

#include "kmap.h"

#include "dissimilarity.h"
#include "warping.h"
#include "utilities.h"
#include "kma_model.h"
#include "checkin.h"

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
                const arma::urowvec& par_opt)
{
    KmaModel model(
        x,
        y,
        n_clust,
        seeds,
        warping_method,
        center_method,
        similarity_method,
        optim_method,
        warping_opt,
        center_opt,
        out_opt(0),
        out_opt(1),
        out_opt(2),
        fence,
        check_total_similarity,
        show_iter,
        comp_original_center,
        par_opt
    );

    return model.execute();
}

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

#ifndef _KMA_MODEL_HPP
#define _KMA_MODEL_HPP

#include <RcppArmadillo.h>

#include "dissimilarity.h"
#include "warping.h"
#include "center_methods.h"
#include "optimizer.h"


using namespace arma;

/// Main class.
/** This class handles loading of the problem and execution of the algorithm.
 */
class KmaModel
{

    const mat x;
    const cube y;
    const uword n_clust;
    const rowvec seeds;

    std::shared_ptr<Dissimilarity> dissim;
    std::shared_ptr<CenterMethod> cen;
    std::shared_ptr<WarpingFunction> warping;
    std::shared_ptr<OptimizerMethod> optimizer;

    const std::string optim_method;
    const rowvec warping_opt;
    const uword n_out;
    const double toll;
    const uword iter_max;
    bool fence;
    bool check_total_similarity;
    bool show_iter;
    bool com_oc;
    uword n_obs;
    uword n_camp;
    uword n_dim;
    const urowvec parallel_opt;

public:
    /// Constructor that load the problem.
    KmaModel(
        const mat& t_x,
        const cube& t_y,
        uword t_n_clust,
        const rowvec& t_seeds,
        std::string warping_method,
        std::string center_method,
        std::string similarity_method,
        std::string t_optim_method,
        const rowvec& t_warping_opt,
        const rowvec& t_center_opt,
        double t_n_out,
        double t_toll,
        uword t_iter_max,
        bool t_fence,
        bool t_check_total_similarity,
        bool t_show_iter,
        bool comp_original_center,
        urowvec par_opt
    );

    /// Method the execute the algorithm.
    Rcpp::List execute();

};

#endif

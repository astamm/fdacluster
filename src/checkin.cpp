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

#include "checkin.h"


#include<RcppArmadillo.h>
#include "dissimilarity.h"
#include "warping.h"
#include "center_methods.h"
#include "optimizer.h"


void checkIn(const arma::mat& x,
             const arma::cube& y,
             const std::string warping_method,
             const std::string center_method,
             const std::string similarity_method,
             const std::string optim_method,
             const util::SharedFactory<WarpingFunction> warfac,
             const util::SharedFactory<Dissimilarity> disfac,
             const util::SharedFactory<CenterMethod> cenfac,
             const util::SharedFactory<OptimizerMethod> optfac,
             const arma::rowvec& warping_opt,
             const arma::rowvec& center_opt,
             const arma::urowvec par_opt,
             bool show_iter)
{

    //
    //  check on x y dimensions
    //
    if(x.n_rows != y.n_rows || x.n_cols != y.n_cols)
        Rcpp::stop(" x dim and y dim must agree");

    //
    //check on available options in factories
    //

    if( warfac.map.count(warping_method)==0)
        Rcpp::stop("warping_method not registered");
    if( disfac.map.count(similarity_method)==0)
        Rcpp::stop("similarity_method not registered");
    if( cenfac.map.count(center_method)==0)
        Rcpp::stop("center_method not registered");
    if( optfac.map.count(optim_method)==0)
        Rcpp::stop("optim_method not registered");

    //
    //  check on warping options
    //

    //
    // check on parallel options
    //

    if(par_opt(1)!=0 & par_opt(1)!=1)
        Rcpp::stop("Parallelization mode could be 0 or 1");

    if(par_opt(1)==1)
        {
            if(center_method!="medoid")
                Rcpp::stop("Parallelization mode 1 is available only with medoid");
        }

    if(show_iter == true)
        Rcpp::Rcout<<"Done"<<endl;
    return;
}

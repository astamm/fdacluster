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

#include <RcppArmadillo.h>
#include "fence.h"
#include "utilities.h"
#include "optimizer.h"


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
    const bool show_iter)
{
    // reading dimension
    arma::uword n_dim = y.n_slices;
    arma::uword n_par = parameters.n_rows;
    arma::urowvec sel_out;
    arma::mat bounds(2,n_par);

    //compute sel_out : observation with outliers
    for(arma::uword i=0; i<n_par; i++)
        {
            double q1=util::quantile(parameters.row(i),0.25);
            double q3=util::quantile(parameters.row(i),0.75);
            bounds(0,i)= q1 - 1.5 * (q3-q1);
            bounds(1,i)= q3 + 1.5 * (q3-q1);

            sel_out= join_horiz( sel_out,
                                 util::which_out(parameters.row(i),bounds(0,i),bounds(1,i))
                               );
        }

    sel_out = unique(sel_out);
    arma::uword it=1;

    // plot FENCE
    if(show_iter==true)
        {
            Rcpp::Rcout<<"IterativeFence: ";
        }


    while(sel_out.n_cols!=0 && it <10)
        {
            if(show_iter==true)
                {
                    Rcpp::Rcout<<" it."<<it<<" ";
                }

            //optimization foe each observation with outliers
            warping->set_bounds(bounds);
            for(uword i=0; i<sel_out.n_cols; i++)
                {
                    uword obs=sel_out[i];
                    uword nt = templates.n_rows;
                    rowvec index_temp(nt);
                    mat parameters_temp(n_par,nt);
                    colvec arg(n_par);
                    mat y_reg = util::approx( x_reg.row(obs), util::observation(y,obs), x_out );
                    for(uword t=0; t<nt; t++ )
                        {
                            mat t_in = templates(span(t),span::all,span::all);
                            if(n_dim >1) t_in =t_in.t();
                            warping_set wset = warping->set_function(x_out, x_out, y_reg, t_in, dissim);
                            auto fun = [&] (const rowvec& arg)
                            {
                                return warping->warp(wset,arg);
                            };

                            index_temp(t) = optimizer->optimize( arg, warping, fun);
                            parameters_temp.col(t) = arg;
                        }// end iterations on templates

                    index(obs) = min( index_temp );
                    labels(obs) = ict( index_min(index_temp));
                    parameters.col(obs)= parameters_temp.col(index_min(index_temp));
                }//end iterations on sel_out


            //recompute sel out
            sel_out.reset();

            for(arma::uword i=0; i<n_par; i++)
                {
                    double q1=util::quantile(parameters.row(i),0.25);
                    double q3=util::quantile(parameters.row(i),0.75);

                    bounds(0,i)= q1 - 1.5 * (q3-q1);
                    bounds(1,i)= q3 + 1.5 * (q3-q1);

                    sel_out= join_horiz( sel_out,
                                         util::which_out(parameters.row(i),bounds(0,i),bounds(1,i))
                                       );
                }
            sel_out = unique(sel_out);
            it++;


        }//fine while

    if(show_iter==true)
        {
            Rcpp::Rcout<<" Done."<<endl;
        }
    return;
}

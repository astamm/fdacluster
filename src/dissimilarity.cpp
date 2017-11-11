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

#include "dissimilarity.h"
#include "utilities.h"

using namespace arma;


grid Dissimilarity::setGrid(const rowvec& xf, const rowvec& xg,
                            const mat& yf, const mat& yg)
{

    grid out;

    out.yf_sim.reset();
    out.yg_sim.reset();
    out.x_sim.reset();

    // dimensions of inputs function
    uword n_dim = yf.n_rows;
    uword nf = xf.size();
    uword ng = xg.size();

    //
    // Na detection and elimination
    //

    mat y_f( n_dim, nf ),y_g( n_dim, ng );
    rowvec x_f(nf), x_g(ng);

    y_f.zeros();
    y_g.zeros();
    x_f.zeros();
    x_g.zeros();


    uword c=0;

    for(uword i=0; i < nf ; i++)
        {
            if( yf.col(i).is_finite() )
                {
                    x_f(c) = xf(i);
                    y_f.col(c) = yf.col(i);
                    c++;
                }
        }

    x_f.resize(c);
    y_f.resize(n_dim,c);

    c=0;
    for(uword i=0; i<ng; i++)
        {
            if( yg.col(i).is_finite() )
                {
                    x_g(c) = xg(i);
                    y_g.col(c) = yg.col(i);
                    c++;
                }
        }


    x_g.resize(c);
    y_g.resize(n_dim,c);

    if( y_f.n_cols == 0 || y_g.n_cols == 0)
        {
            Rcpp::warning("Una delle due funzioni è completamente NA.");
            return out;
        }

////////////////////// x_sim e approssimazioni /////////////////////////////////////////

    double x_min = std::max( x_f.min(), x_g.min() );
    double x_max = std::min( x_f.max(), x_g.max() );

    double pf(0),pg(0);
    for (uword i=0; i<x_f.size(); i++)
        {
            if(x_f(i)>=x_min && x_f(i)<=x_max)
                pf++;
        }
    for (uword i=0; i<x_g.size(); i++)
        {
            if(x_g(i)>=x_min && x_g(i)<=x_max)
                pg++;
        }

    double p = std::max(pg,pf);

    if( p<=1 )
        {
            Rcpp::warning("Disimilarity: no open intervals in common. Reduction of warping bounds suggested.");
            return out;
        }


    rowvec xs(p);
    xs.zeros();

    xs(0) = x_min;
    for (size_t i = 1; i<p; i++)
        {
            xs(i) =  xs(i-1) + (x_max-x_min) / (p-1);
        }

    out.yf_sim.set_size(n_dim,p);
    out.yg_sim.set_size(n_dim,p);

    // set attributes of dissimilarity
    out.yf_sim = util::approx( x_f, y_f, xs );
    out.yg_sim = util::approx( x_g, y_g, xs );
    out.x_sim = xs;
    return out;
}


double Pearson::compute(const rowvec& xf, const rowvec& xg,
                        const mat& yf, const mat& yg)
{
    grid gr = setGrid(xf,xg,yf,yg);

    if(gr.yf_sim.is_empty())
        return 10000000;

    uword dim = gr.yf_sim.n_rows;
    double res(0);

    for(uword k=0; k < dim; k++ )
        {
            res += norm_dot( gr.yf_sim.row(k), gr.yg_sim.row(k) );
        }
    return -res/dim;
}



double L2::compute(const rowvec& xf, const rowvec& xg,
                   const mat& yf, const mat& yg)
{
    grid gr=setGrid(xf,xg,yf,yg);
    if(gr.yf_sim.is_empty())
        return 10000000;

    uword dim = gr.yf_sim.n_rows;
    uword len = gr.x_sim.size();

    double res(0);

    rowvec d = gr.x_sim.cols(1,len-1) - gr.x_sim.cols(0,len-2);
    double D = gr.x_sim(len-1)-gr.x_sim(0);

    for(uword k=0; k < dim; k++ )
        {
            rowvec diff = sqrt(d) % ( gr.yf_sim.row(k).cols(1,len-1)- gr.yg_sim.row(k).cols(1,len-1) );
            res += dot(diff,diff)/(D*dim);
        }

    return sqrt(res);

}


double L2w::compute(const rowvec& xf, const rowvec& xg,
                   const mat& yf, const mat& yg)
{
  grid gr=setGrid(xf,xg,yf,yg);
  if(gr.yf_sim.is_empty())
    return 10000000;

  uword dim = gr.yf_sim.n_rows;
  uword len = gr.x_sim.size();

  double res(0);

  rowvec d = gr.x_sim.cols(1,len-1) - gr.x_sim.cols(0,len-2);
  double D = gr.x_sim(len-1)-gr.x_sim(0);

  rowvec w(len-1);
  for(uword i=0;i < (len-1); i++)
    w(i)=1/(i+1);

  for(uword k=0; k < dim; k++ )
  {
    rowvec diff =  sqrt(w % d) % ( gr.yf_sim.row(k).cols(1,len-1)- gr.yg_sim.row(k).cols(1,len-1) );
    res += dot(diff, diff)/(D*dim);
  }

  return sqrt(res);

}


double L2first::compute(const rowvec& xf, const rowvec& xg,
                    const mat& yf, const mat& yg)
{
  grid gr=setGrid(xf,xg,yf,yg);
  if(gr.yf_sim.is_empty())
    return 10000000;

  uword dim = gr.yf_sim.n_rows;
  uword len = gr.x_sim.size();

  double res(0);

  rowvec d = gr.x_sim.cols(1,len-1) - gr.x_sim.cols(0,len-2);
  double D = gr.x_sim(len-1)-gr.x_sim(0);

  rowvec w= Rcpp::rep(0.001,len-1);
  w(0)=1;

  for(uword k=0; k < dim; k++ )
  {
    rowvec diff =  sqrt(w % d) % ( gr.yf_sim.row(k).cols(1,len-1)- gr.yg_sim.row(k).cols(1,len-1) );
    res += dot(diff, diff)/(D*dim);
  }

  return sqrt(res);

}

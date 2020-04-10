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

#include "affineWarpingClass.h"

arma::mat AffineWarpingFunction::ApplyWarping(const arma::mat &x, const arma::mat &par)
{
    arma::mat out(x.n_rows, x.n_cols);

    for (unsigned int i = 0;i < x.n_rows;++i)
        out.row(i) = par(0, i) * x.row(i) + par(1, i);

    return out;
}

void AffineWarpingFunction::SetParameterBounds(const arma::rowvec &war_opt,
                                               const arma::mat &x)
{
    double workScalar = arma::as_scalar(arma::min(arma::max(x, 1) - arma::min(x, 1)));
    double dl = war_opt(0), sl = war_opt(1);
    m_ParameterLowerBounds = { 1 - dl, -sl * workScalar};
    m_ParameterUpperBounds = { 1 + dl,  sl * workScalar};
}

arma::mat AffineWarpingFunction::GetFinalWarping(const arma::cube &parameters_vec,
                                                 const arma::urowvec &labels,
                                                 const arma::urowvec &ict)
{
    unsigned int np = parameters_vec.n_rows;
    unsigned int n_obs = parameters_vec.n_cols;
    unsigned int n_iter = parameters_vec.n_slices;
    arma::mat out(np,n_obs);
    out.row(0).ones();
    arma::rowvec a;
    arma::rowvec b;

    for (unsigned int i = 0;i < n_iter;++i)
    {
        a = parameters_vec(arma::span(0), arma::span::all, arma::span(i));
        b = parameters_vec(arma::span(1), arma::span::all, arma::span(i));
        out.row(0) = out.row(0) % a;
        out.row(1) = out.row(1) % a + b;
    }

    arma::urowvec sel;
    arma::colvec m;

    for (unsigned int k = 0;k < ict.size();++k)
    {
        sel = find(labels == ict(k)).t();

        //compute means
        m = sum(out.cols(sel), 1) / sel.size();

        for (unsigned int i = 0;i < sel.size();++i)
        {
            out(0, sel(i)) = out(0, sel(i)) / m(0);
            out(1, sel(i)) = (out(1, sel(i)) - m(1)) / m(0);
        }
    }

    return out;
}

void AffineWarpingFunction::Normalize(arma::mat &par,
                                      const arma::urowvec &ict,
                                      const arma::urowvec &labels)
{
    arma::urowvec sel;
    arma::colvec par_mean;

    for (unsigned int i = 0;i < ict.size();++i)
    {
        sel = find(labels == ict(i)).t();

        // calcolo medie cluster
        par_mean = mean(par.cols(sel), 1);

        // aggiorno shift e dilation
        for (unsigned int j = 0;j < sel.size();++j)
        {
            // normalized dilation
            par(0, sel(j)) = par(0, sel(j)) / par_mean(0);
            // normalized shift
            par(1, sel(j)) = -par(0, sel(j)) * par_mean(1) / par_mean(0) + par(1, sel(j));
        }
    }
}

double AffineWarpingFunction::GetDissimilarityAfterWarping(const WarpingSet &w_set,
                                                           const arma::colvec &arg)
{
    return m_DissimilarityPointer->GetDistance(arg(0) * w_set.xf + arg(1), w_set.xg, w_set.yf, w_set.yg);
}

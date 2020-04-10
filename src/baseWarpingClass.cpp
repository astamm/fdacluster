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

#include "baseWarpingClass.h"

WarpingSet BaseWarpingFunction::SetInputData(const arma::rowvec &x_f,
                                         const arma::rowvec &x_g,
                                         const arma::mat &y_f,
                                         const arma::mat &y_g,
                                         const std::shared_ptr<Dissimilarity> &d)
{
    WarpingSet out;

    out.xf = x_f;
    out.xg = x_g;
    out.yf = y_f;
    out.yg = y_g;

    m_DissimilarityPointer = d;

    return out;
}

//-------------------
// Dilation warping
//-------------------

void DilationFunction::normalize(mat& par,const urowvec& ict,const urowvec& labels)
{

    for(size_t i=0; i<ict.size(); i++)
        {
            urowvec sel =find(labels==ict(i)).t();

            //calcolo medie cluster
            colvec par_mean= mean(par.cols(sel),1);

            // aggiorno shift e dilation
            for(size_t j=0; j<sel.size(); j++)
                {
                    par(0,sel(j))=  par(0,sel(j))/par_mean(0);
                }// for su sel

        }// for su ict
}
void ShiftFunction::normalize(mat& par,const urowvec& ict,const urowvec& labels)
{

    for(size_t i=0; i<ict.size(); i++)
        {
            urowvec sel =find(labels==ict(i)).t();

            //calcolo medie cluster
            colvec par_mean= mean(par.cols(sel),1);

            // aggiorno shift e dilation
            for(size_t j=0; j<sel.size(); j++)
                {

                    par(0,sel(j))=  par(0,sel(j)) -par_mean(0);

                }// for su sel
        }// for su ict
}

//
// Final Warping
//

mat ShiftFunction::final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict)
{
    uword n_obs = parameters_vec.n_cols;
    uword n_iter = parameters_vec.n_slices;
    mat out(1,n_obs);
    out.row(0).zeros();
    for(size_t i=0; i<n_iter; i++)
        {
            rowvec a = parameters_vec(span(0),span::all,span(i));
            out.row(0) = out.row(0) + a;
        }
    for(uword k = 0; k < ict.size(); k++)
        {
            urowvec sel = find(labels == ict(k)).t();
            colvec m = mean(out.cols(sel),1);
            for(uword i =0; i< sel.size(); i++)
                {
                    out(0,sel(i)) = (out(0,sel(i)) - m(0));
                }
        }
    return out;
}
mat DilationFunction::final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict)
{
    uword n_obs = parameters_vec.n_cols;
    uword n_iter = parameters_vec.n_slices;
    mat out(1,n_obs);
    out.row(0).ones();
    for(size_t i=0; i<n_iter; i++)
        {
            rowvec a = parameters_vec(span(0),span::all,span(i));
            out.row(0) = out.row(0) % a ;
        }
    for(uword k = 0; k < ict.size(); k++)
        {
            urowvec sel = find(labels == ict(k)).t();
            //compute means
            colvec m = mean(out.cols(sel),1);
            for(size_t i =0; i< sel.size(); i++)
                {
                    out(0,sel(i)) = out(0,sel(i))/m(0);
                }
        }
    return out;
}
mat NoAlignmentFunction::final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict)
{
    mat out(0,labels.n_cols);
    return out;

}


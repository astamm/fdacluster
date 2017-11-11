// Copyright (C) 2017 Alessandro Zito  (zito.ales@gmail.com)
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
#ifdef _OPENMP
 #include <omp.h>
#endif


#include "newcenters.h"
#include "warping.h"
#include "center_methods.h"

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
){
  // switch to choose how to parallalelize
  // case 0 trivial  each thread one cluster
  //case 1 each cluster all the threads (available only with medoid)

  switch(par_opt(1))
  {
  case 0:

    #ifdef _OPENMP
      #pragma omp parallel for num_threads(par_opt(0))
    #endif
    for(uword i=0; i< ict.size(); i++)
    {
      urowvec sel = find(labels == ict(i)).t();
      center a = cen->computeCenter( x_reg.rows(sel), util::observations(y,sel), dissim, x_out );
      templates.tube(span(i),span::all) = a.y_center.t();

// #ifdef _OPENMP
// #pragma omp critical
// {
//   if(show_iter==true)
//     cout<<"The thread num. "<<omp_get_thread_num()<<" has computed template num. "<<i<<endl;
// }
// #endif
    }
    break;

  case 1:

    for(uword i=0; i< ict.size(); i++)
        {
            urowvec sel = find(labels == ict(i)).t();
            center a = cen->computeParallelCenter( x_reg.rows(sel), util::observations(y,sel), dissim, x_out, par_opt(0));
            templates.tube(span(i),span::all) = a.y_center.t();
            if(show_iter==true)
                cout<<"Template num. "<<i<<" computed"<<endl;
        }
    break;
  }

}// fine new centers

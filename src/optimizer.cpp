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
#include <memory>

#include "optimizer.h"
#include "dlib/optimization/find_optimal_parameters.h"

typedef dlib::matrix<double,0,1> argument;

double Bobyqa::optimize(colvec& arg, std::shared_ptr<WarpingFunction>& pfunc, std::function<double(colvec)> fun)
{
    //LAMBDA FUNCTION
    auto fun2 = [&fun] (const argument& argt)
    {
        // here I have to convert the argument input to an input for fun
        colvec argt_fun(argt.nr());

        for(uword i=0; i<argt.nr(); i++)
            argt_fun(i)=argt(i);
        return fun(argt_fun);
    };

    uword dp = pfunc->n_pars();
    rowvec s= (pfunc->get_lower_bound() + pfunc->get_upper_bound())/2;
    rowvec d= ( pfunc->get_upper_bound() - pfunc->get_lower_bound());

    double eps = 0.00000001;

    argument starting_point(dp);
    argument lbound(dp);
    argument ubound(dp);

    for(uword i =0; i < dp; i++)
        {
            starting_point(i)=s(i);
            lbound(i)=(pfunc->get_lower_bound())(i);
            ubound(i)=(pfunc->get_upper_bound())(i);
        }

    //gestione caso no aligment numero di parametri da stimare 0
    if(arg.size()==0)
        return fun2(starting_point);


    double radius = min(d)/2-0.0000001;

    find_optimal_parameters (
        radius,
        eps,
        100,
        starting_point,
        lbound,
        ubound,
        fun2
    );

    for(uword i =0; i < dp; i++)
        {
            arg(i) = starting_point(i);
        }

    return fun2(starting_point);

}



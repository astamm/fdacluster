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

#ifndef MY_OPTIMIZER_HPP
#define MY_OPTIMIZER_HPP

#include <RcppArmadillo.h>
#include "warping.h"

//
//  optimizer
//

/// Optimizer Base class
class OptimizerMethod
{
public:
    ///optimer member function
    /**
     * @param[arg] warping's warping to be optimized returned by reference;
     * @param[pfunc] pointer to warping object;
     * @param[fun] functor(colvec) to be optimized, input customizable.
     *
     * return optimal dissimilarity.
     */
    virtual double optimize(colvec& arg, std::shared_ptr<WarpingFunction>& pfunc, std::function<double(colvec)> fun)=0;
};


class Bobyqa: public OptimizerMethod
{
public:
    virtual double optimize(colvec& arg, std::shared_ptr<WarpingFunction>& pfunc, std::function<double(colvec)> fun);
};

#endif

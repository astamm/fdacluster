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

#ifndef DISSIMILARITY_HPP_
#define DISSIMILARITY_HPP_

#include<RcppArmadillo.h>

using namespace arma;

/// Element with approximated abiscissa and approximated functions to compute dissimilarity.
struct grid
{
    mat yf_sim;
    mat yg_sim;
    rowvec x_sim;
};


/// Base class for all the available dissimilarity.
class Dissimilarity
{
public:

    Dissimilarity() {};
    /// compute dissimilarity method different for each derived class
    /**
     * @param[xf] abscissa of the f function;
     * @param[xg] abscissa of the g function;
     * @param[yf] values of the f function;
     * @param[yg] values of the g function;
     *
     * @return dissimilarity between f and g functions.
     */
    virtual double compute(const rowvec& xf, const rowvec& xg,const mat& yf, const mat& yg)=0;

    /// compute common grid and approximations to compute dissimilarity
    /**
    * @param[xf] abscissa of the f function;
    * @param[xg] abscissa of the g function;
    * @param[yf] values of the f function;
    * @param[yg] values of the g function;
    *
    * @return grid element with approximated functions on a common grid.
    */
    grid setGrid(const rowvec& xf, const rowvec& xg, const mat& yf, const mat& yg);
};


/// Pearson similarity
/**
  *The function is returned negative to be a dissimilarity and
  * fit the minimization of the algorithm.
  */
class Pearson final: public Dissimilarity
{
public:

    Pearson():Dissimilarity() {};

    virtual double compute(const rowvec& xf, const rowvec& xg,
                           const mat& yf, const mat& yg);
};

/// L2 Distance
class L2 final: public Dissimilarity
{
public:
    L2():Dissimilarity() {};
    virtual double compute(const rowvec& xf, const rowvec& xg,
                           const mat& yf, const mat& yg);
};

/// L2w Distance
class L2w final: public Dissimilarity
{
public:
  L2w():Dissimilarity() {};
  virtual double compute(const rowvec& xf, const rowvec& xg,
                         const mat& yf, const mat& yg);
};


/// L2w Distance
class L2first final: public Dissimilarity
{
public:
  L2first():Dissimilarity() {};
  virtual double compute(const rowvec& xf, const rowvec& xg,
                         const mat& yf, const mat& yg);
};

#endif

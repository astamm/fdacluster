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

#include <RcppArmadillo.h>

/// Element with approximated abiscissa and approximated functions to compute dissimilarity.
struct FunctionPair
{
    arma::rowvec Grid;
    arma::mat Values1;
    arma::mat Values2;
};

/// Base class for all the available dissimilarity.
class Dissimilarity
{
public:
    /// compute dissimilarity method different for each derived class
    /**
     * @param[grid1] evaluation grid of the first function;
     * @param[grid2] evaluation grid of the second function;
     * @param[values1] values of the first function;
     * @param[values2] values of the second function;
     *
     * @return dissimilarity between the two input functions.
     */
    virtual double GetDistance(const arma::rowvec& grid1,
                               const arma::rowvec& grid2,
                               const arma::mat& values1,
                               const arma::mat& values2) = 0;

    /// compute common grid and approximations to compute dissimilarity
    /**
     * @param[grid1] evaluation grid of the first function;
     * @param[grid2] evaluation grid of the second function;
     * @param[values1] values of the first function;
     * @param[values2] values of the second function;
    *
    * @return FunctionPair object containing the two functions evaluated on their common grid.
    */
    FunctionPair GetComparableFunctions(const arma::rowvec& grid1,
                                        const arma::rowvec& grid2,
                                        const arma::mat& values1,
                                        const arma::mat& values2);

protected:
    Dissimilarity() {}
    virtual ~Dissimilarity() {}
};

/// Pearson similarity
/**
  *The function is returned negative to be a dissimilarity and
  * fit the minimization of the algorithm.
  */
class Pearson : public Dissimilarity
{
public:
    virtual double GetDistance(const arma::rowvec& grid1,
                               const arma::rowvec& grid2,
                               const arma::mat& values1,
                               const arma::mat& values2);
};

/// L2 Distance
class L2 : public Dissimilarity
{
public:
    virtual double GetDistance(const arma::rowvec& grid1,
                               const arma::rowvec& grid2,
                               const arma::mat& values1,
                               const arma::mat& values2);
};

/// L2w Distance
class L2w : public Dissimilarity
{
public:
    virtual double GetDistance(const arma::rowvec& grid1,
                               const arma::rowvec& grid2,
                               const arma::mat& values1,
                               const arma::mat& values2);
};


/// L2w Distance
class L2first : public Dissimilarity
{
public:
    virtual double GetDistance(const arma::rowvec& grid1,
                               const arma::rowvec& grid2,
                               const arma::mat& values1,
                               const arma::mat& values2);
};

#endif

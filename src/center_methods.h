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

#ifndef CENTER_METHODS_HPP_
#define CENTER_METHODS_HPP_

#include "dissimilarity.h"

#include <RcppArmadillo.h>
#include <memory>

/// Center element returned by computeCenter and computeParallelCenter methods
struct CenterObject
{
    arma::rowvec Grid;
    arma::mat Values;
    arma::rowvec Distances;
};


/// Base class for all the center methods available.
class CenterMethod
{
public:
    /// Compute center method.
    /**
     *  @param[inputGrid] Input grid on which observed functions are evaluated;
     *  @param[inputValues] Input function values on input grid;
     *  @param[distanceObject] Shared pointer to the base class Dissimilarity;
     *  @param[outputGrid] Output grid on which the center will be evaluated.
     *
     *  @return A center object.
     */
    virtual CenterObject GetCenter(const arma::mat& inputGrid,
                                   const arma::cube& inputValues,
                                   std::shared_ptr<Dissimilarity>& distanceObject,
                                   const arma::rowvec& outputGrid) = 0;

    /// Compute center method in parallel (used if type of parallelization is 1).
    /**
     *  @param[inputGrid] Input grid on which observed functions are evaluated;
     *  @param[inputValues] Input function values on input grid;
     *  @param[distanceObject] Shared pointer to the base class Dissimilarity;
     *  @param[outputGrid] Output grid on which the center will be evaluated;
     *  @param[nbThreads] Number of threads to use during the computation.
     *
     *  @return A center object.
     */
    virtual CenterObject GetCenterParallel(const arma::mat& inputGrid,
                                           const arma::cube& inputValues,
                                           std::shared_ptr<Dissimilarity>& distanceObject,
                                           const arma::rowvec& outputGrid,
                                           unsigned int nbThreads)
    {
        CenterObject outputCenter;
        return outputCenter;
    }

    void SetSpanValue(const double num) {m_SpanValue = num;}
    double GetSpanValue() {return m_SpanValue;}

    void SetDeltaValue(const double num) {m_DeltaValue = num;};
    double GetDeltaValue() {return m_DeltaValue;}

protected:
    CenterMethod()
    {
        m_SpanValue = 0.0;
        m_DeltaValue = 0.0;
    }

    virtual ~CenterMethod() {};

private:
    double m_SpanValue, m_DeltaValue;
};

/// Medoid center method finds the real medoid center.
class Medoid : public CenterMethod
{
public:
    virtual CenterObject GetCenter(const arma::mat& inputGrid,
                                   const arma::cube& inputValues,
                                   std::shared_ptr<Dissimilarity>& distanceObject,
                                   const arma::rowvec& outputGrid);

    virtual CenterObject GetCenterParallel(const arma::mat& inputGrid,
                                           const arma::cube& inputValues,
                                           std::shared_ptr<Dissimilarity>& distanceObject,
                                           const arma::rowvec& outputGrid,
                                           unsigned int nbThreads);
};


/// PseudoMedoid center method finds the real medoid of each component and create a center.
class PseudoMedoid : public CenterMethod
{
public:
    virtual CenterObject GetCenter(const arma::mat& inputGrid,
                                   const arma::cube& inputValues,
                                   std::shared_ptr<Dissimilarity>& distanceObject,
                                   const arma::rowvec& outputGrid);
};

/// Mean center method compute an approximation by lowess.
class Mean : public CenterMethod
{
public:
    virtual CenterObject GetCenter(const arma::mat& inputGrid,
                                   const arma::cube& inputValues,
                                   std::shared_ptr<Dissimilarity>& distanceObject,
                                   const arma::rowvec& outputGrid);
};

/// Median center method.
class Median : public CenterMethod
{
public:
    virtual CenterObject GetCenter(const arma::mat& inputGrid,
                                   const arma::cube& inputValues,
                                   std::shared_ptr<Dissimilarity>& distanceObject,
                                   const arma::rowvec& outputGrid);
};

#endif

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
#include <squad.h>

#include "dissimilarity.h"
#include "utilities.h"

FunctionPair Dissimilarity::GetComparableFunctions(const arma::rowvec& grid1,
                                                   const arma::rowvec& grid2,
                                                   const arma::mat& values1,
                                                   const arma::mat& values2)
{
    FunctionPair outputPair;
    outputPair.Grid.reset();
    outputPair.Values1.reset();
    outputPair.Values2.reset();

    // dimensions of inputs function
    unsigned int nDim = values1.n_rows;

    if (values2.n_rows != nDim)
        Rcpp::stop("Function domains have not the same dimension.");

    unsigned int nPts1 = values1.n_cols;

    if (grid1.size() != nPts1)
        Rcpp::stop("Function 1 is not properly evaluated on its grid.");

    unsigned int nPts2 = values2.n_cols;

    if (grid2.size() != nPts2)
        Rcpp::stop("Function 2 is not properly evaluated on its grid.");

    //
    // Missing value detection and elimination
    //

    arma::mat cleanValues1(nDim, nPts1, arma::fill::zeros), cleanValues2(nDim, nPts2, arma::fill::zeros);
    arma::rowvec cleanGrid1(nPts1, arma::fill::zeros), cleanGrid2(nPts2, arma::fill::zeros);

    unsigned int c = 0;
    for (unsigned int i = 0;i < nPts1;++i)
    {
        if (arma::is_finite(grid1(i)) && arma::is_finite(values1.col(i)))
        {
            cleanGrid1(c) = grid1(i);
            cleanValues1.col(c) = values1.col(i);
            ++c;
        }
    }

    cleanGrid1.resize(c);
    cleanValues1.resize(nDim, c);

    c = 0;
    for (unsigned int i = 0;i < nPts2;++i)
    {
        if (arma::is_finite(grid2(i)) && arma::is_finite(values2.col(i)))
        {
            cleanGrid2(c) = grid2(i);
            cleanValues2.col(c) = values2.col(i);
            ++c;
        }
    }

    cleanGrid2.resize(c);
    cleanValues2.resize(nDim, c);

    if (cleanValues1.n_cols == 0 || cleanValues2.n_cols == 0)
    {
        Rcpp::warning("Dissimilarity: at least one function only contains missing values.");
        return outputPair;
    }

    // Compute smallest common grid

    double xMin = std::max(cleanGrid1.min(), cleanGrid2.min());
    double xMax = std::min(cleanGrid1.max(), cleanGrid2.max());

    if (xMin >= xMax)
    {
        Rcpp::warning("Disimilarity: domain intersection is empty.");
        return outputPair;
    }

    arma::rowvec xCommon = arma::linspace<arma::rowvec>(xMin, xMax, (nPts1 + nPts2) / 2);

    outputPair.Grid = xCommon;
    outputPair.Values1 = util::approx(cleanGrid1, cleanValues1, xCommon);
    outputPair.Values2 = util::approx(cleanGrid2, cleanValues2, xCommon);

    return outputPair;
}

FunctionPair Dissimilarity::GetComparableQuaternionFunctions(const arma::rowvec& grid1,
                                                             const arma::rowvec& grid2,
                                                             const arma::mat& values1,
                                                             const arma::mat& values2)
{
    FunctionPair outputPair;
    outputPair.Grid.reset();
    outputPair.Values1.reset();
    outputPair.Values2.reset();

    // dimensions of inputs function
    unsigned int nDim = values1.n_rows;

    if (values2.n_rows != nDim)
        Rcpp::stop("Function domains have not the same dimension.");

    unsigned int nPts1 = values1.n_cols;

    if (grid1.size() != nPts1)
        Rcpp::stop("Function 1 is not properly evaluated on its grid.");

    unsigned int nPts2 = values2.n_cols;

    if (grid2.size() != nPts2)
        Rcpp::stop("Function 2 is not properly evaluated on its grid.");

    //
    // Missing value detection and elimination
    //

    arma::mat cleanValues1(nDim, nPts1, arma::fill::zeros), cleanValues2(nDim, nPts2, arma::fill::zeros);
    arma::rowvec cleanGrid1(nPts1, arma::fill::zeros), cleanGrid2(nPts2, arma::fill::zeros);

    unsigned int c = 0;
    for (unsigned int i = 0;i < nPts1;++i)
    {
        if (arma::is_finite(grid1(i)) && arma::is_finite(values1.col(i)))
        {
            cleanGrid1(c) = grid1(i);
            cleanValues1.col(c) = values1.col(i);
            ++c;
        }
    }

    cleanGrid1.resize(c);
    cleanValues1.resize(nDim, c);

    c = 0;
    for (unsigned int i = 0;i < nPts2;++i)
    {
        if (arma::is_finite(grid2(i)) && arma::is_finite(values2.col(i)))
        {
            cleanGrid2(c) = grid2(i);
            cleanValues2.col(c) = values2.col(i);
            ++c;
        }
    }

    cleanGrid2.resize(c);
    cleanValues2.resize(nDim, c);

    if (cleanValues1.n_cols == 0 || cleanValues2.n_cols == 0)
    {
        Rcpp::warning("Dissimilarity: at least one function only contains missing values.");
        return outputPair;
    }

    // Compute smallest common grid

    double xMin = std::max(cleanGrid1.min(), cleanGrid2.min());
    double xMax = std::min(cleanGrid1.max(), cleanGrid2.max());

    if (xMin >= xMax)
    {
        Rcpp::warning("Disimilarity: domain intersection is empty.");
        return outputPair;
    }

    unsigned int nPts = std::min(cleanGrid1.size(), cleanGrid2.size());

    outputPair.Grid = arma::linspace<arma::rowvec>(xMin, xMax, nPts);

    outputPair.Values1 = Rcpp::as<arma::mat>(squad::RegularizeGrid(
        Rcpp::wrap(cleanGrid1),
        Rcpp::wrap(cleanValues1),
        xMin, xMax, nPts
    ));

    outputPair.Values2 = Rcpp::as<arma::mat>(squad::RegularizeGrid(
        Rcpp::wrap(cleanGrid2),
        Rcpp::wrap(cleanValues2),
        xMin, xMax, nPts
    ));

    return outputPair;
}

double Pearson::GetDistance(const arma::rowvec& grid1,
                            const arma::rowvec& grid2,
                            const arma::mat& values1,
                            const arma::mat& values2)
{
    FunctionPair pair = this->GetComparableFunctions(grid1, grid2, values1, values2);

    if (pair.Grid.is_empty())
        return 1e7;

    unsigned int nDim = pair.Values1.n_rows;
    double pearsonCorrelation = 0.0;

    for (unsigned int k = 0;k < nDim;++k)
        pearsonCorrelation += arma::norm_dot(pair.Values1.row(k), pair.Values2.row(k));

    return -pearsonCorrelation / (double)nDim;
}

double L2::GetDistance(const arma::rowvec& grid1,
                       const arma::rowvec& grid2,
                       const arma::mat& values1,
                       const arma::mat& values2)
{
    FunctionPair pair = this->GetComparableFunctions(grid1, grid2, values1, values2);

    if (pair.Grid.is_empty())
        return std::numeric_limits<double>::max();

    unsigned int nDim = pair.Values1.n_rows;
    unsigned int nPts = pair.Grid.size();

    if (nPts <= 1.0)
        return std::numeric_limits<double>::max();

    arma::rowvec diffVector = pair.Grid.cols(1, nPts - 1) - pair.Grid.cols(0, nPts - 2);
    double rangeValue = pair.Grid(nPts - 1) - pair.Grid(0);
    arma::rowvec workVector;
    double squaredDistanceValue = 0.0;

    for (unsigned int k = 0;k < nDim;++k)
    {
        workVector = arma::sqrt(diffVector) % (pair.Values1.row(k).cols(1, nPts - 1) - pair.Values2.row(k).cols(1, nPts - 1));
        squaredDistanceValue += arma::dot(workVector, workVector);
    }

    squaredDistanceValue /= (rangeValue * (double)nDim);

    return std::sqrt(squaredDistanceValue);
}

double UnitQuaternionL2::GetDistance(const arma::rowvec& grid1,
                                     const arma::rowvec& grid2,
                                     const arma::mat& values1,
                                     const arma::mat& values2)
{
    FunctionPair pair = this->GetComparableQuaternionFunctions(
        grid1, grid2, values1, values2
    );

    if (pair.Grid.is_empty())
        return std::numeric_limits<double>::max();

    unsigned int nPts = pair.Grid.size();

    if (nPts <= 1.0)
        return std::numeric_limits<double>::max();

    double squaredDistanceValue = 0.0;

    for (unsigned int j = 0;j < nPts - 1;++j)
    {
        double tmpDistance = squad::GeodesicQuaternionDistance(
            Rcpp::wrap(pair.Values1),
            Rcpp::wrap(pair.Values2),
            j + 1, j + 1
        );

        squaredDistanceValue += tmpDistance * tmpDistance;
    }

    squaredDistanceValue /= (nPts - 1.0);

    return std::sqrt(squaredDistanceValue);
}

double L2w::GetDistance(const arma::rowvec& grid1,
                        const arma::rowvec& grid2,
                        const arma::mat& values1,
                        const arma::mat& values2)
{
    FunctionPair pair = this->GetComparableFunctions(grid1, grid2, values1, values2);

    if (pair.Grid.is_empty())
        return 1e7;

    unsigned int nDim = pair.Values1.n_rows;
    unsigned int nPts = pair.Grid.size();

    if (nPts <= 1.0)
        return 1e7;

    double squaredDistanceValue = 0.0;

    arma::rowvec diffVector = pair.Grid.cols(1, nPts - 1) - pair.Grid.cols(0, nPts - 2);
    double rangeValue = pair.Grid(nPts - 1) - pair.Grid(0);
    arma::rowvec workVector;

    arma::rowvec weightVector(nPts - 1);
    for (unsigned int i = 0;i < nPts - 1;++i)
        weightVector(i) = 1.0 / ((double)i + 1.0);

    for (unsigned int k = 0;k < nDim;++k)
    {
        workVector = arma::sqrt(weightVector % diffVector) % (pair.Values1.row(k).cols(1, nPts - 1) - pair.Values2.row(k).cols(1, nPts - 1));
        squaredDistanceValue += arma::dot(workVector, workVector) / (rangeValue * (double)nDim);
    }

    return std::sqrt(squaredDistanceValue);
}

double L2first::GetDistance(const arma::rowvec& grid1,
                            const arma::rowvec& grid2,
                            const arma::mat& values1,
                            const arma::mat& values2)
{
    FunctionPair pair = this->GetComparableFunctions(grid1, grid2, values1, values2);

    if (pair.Grid.is_empty())
        return 1e7;

    unsigned int nDim = pair.Values1.n_rows;
    unsigned int nPts = pair.Grid.size();

    if (nPts <= 1.0)
        return 1e7;

    double squaredDistanceValue = 0.0;

    arma::rowvec diffVector = pair.Grid.cols(1, nPts - 1) - pair.Grid.cols(0, nPts - 2);
    double rangeValue = pair.Grid(nPts - 1) - pair.Grid(0);
    arma::rowvec workVector;

    arma::rowvec weightVector(nPts - 1);
    weightVector.fill(1e-3);
    weightVector(0) = 1.0;

    for (unsigned int k = 0;k < nDim;++k)
    {
        workVector = arma::sqrt(weightVector % diffVector) % (pair.Values1.row(k).cols(1, nPts - 1) - pair.Values2.row(k).cols(1, nPts - 1));
        squaredDistanceValue += arma::dot(workVector, workVector) / (rangeValue * (double)nDim);
    }

    return std::sqrt(squaredDistanceValue);
}

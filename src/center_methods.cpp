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
#ifdef _OPENMP
#include <omp.h>
#endif

#include "center_methods.h"
#include "low.h"
#include "utilities.h"

CenterObject Medoid::GetCenter(const arma::mat& inputGrid,
                               const arma::cube& inputValues,
                               std::shared_ptr<Dissimilarity>& distanceObject,
                               const arma::rowvec& outputGrid)
{
    CenterObject outputCenter;

    unsigned int nOut = outputGrid.size();
    unsigned int nDim = inputValues.n_slices;
    unsigned int nObs = inputValues.n_rows;
    unsigned int nPts = inputValues.n_cols;

    if (inputGrid.n_cols != nPts)
        Rcpp::stop("The number of columns in x should match the second dimension of y.");

    if (inputGrid.n_rows != nObs)
        Rcpp::stop("The number of rows in x should match the first dimension of y.");

    arma::mat D(nObs, nObs, arma::fill::zeros);
    arma::mat workMatrix1, workMatrix2;

    for (unsigned int i = 0;i < nObs;++i)
    {
        workMatrix1 = util::GetObservation(inputValues, i);

        for (unsigned int j = i+1;j < nObs;++j)
        {
            workMatrix2 = util::GetObservation(inputValues, j);

            D(i, j) = distanceObject->GetDistance(inputGrid.row(i), inputGrid.row(j), workMatrix1, workMatrix2);
            D(j, i) = D(i, j);
        }
    }

    arma::colvec distVec = arma::sum(D, 1);
    unsigned int medoidIndex = arma::index_min(distVec);
    arma::mat medoidValue = util::GetObservation(inputValues, medoidIndex);

    outputCenter.Grid = outputGrid;
    outputCenter.Values = util::approx(inputGrid.row(medoidIndex), medoidValue, outputGrid);
    outputCenter.Distances = D.row(medoidIndex);

    return outputCenter;
}

CenterObject Mean::GetCenter(const arma::mat& inputGrid,
                             const arma::cube& inputValues,
                             std::shared_ptr<Dissimilarity>& distanceObject,
                             const arma::rowvec& outputGrid)
{
    CenterObject outputCenter;

    unsigned int nOut = outputGrid.size();
    unsigned int nDim = inputValues.n_slices;
    unsigned int nObs = inputValues.n_rows;
    unsigned int nPts = inputValues.n_cols;

    if (inputGrid.n_cols != nPts)
        Rcpp::stop("The number of columns in x should match the second dimension of y.");

    if (inputGrid.n_rows != nObs)
        Rcpp::stop("The number of rows in x should match the first dimension of y.");

    arma::mat yOut(nDim, nOut);
    std::vector<double> yIn, xIn;
    std::vector< std::pair<double,double> > zipped;
    std::vector<double> ys;
    arma::rowvec x1, y1;

    for (unsigned int k = 0;k < nDim;++k)
    {
        xIn.clear();
        yIn.clear();
        zipped.clear();
        ys.clear();

        for (unsigned int i = 0;i < nObs;++i)
        {
            for (unsigned int j = 0;j < nPts;++j)
            {
                if (arma::is_finite(inputGrid(i, j)) && arma::is_finite(inputValues(i, j, k)))
                {
                    xIn.push_back(inputGrid(i, j));
                    yIn.push_back(inputValues(i, j, k));
                }
            }
        }

        util::zip(xIn, yIn, zipped);

        // Sort the vector of pairs
        std::sort(std::begin(zipped), std::end(zipped),
                  [&](const std::pair<double,double>& a, const std::pair<double,double>& b)
                  {
                      return a.first < b.first;
                  });

        util::unzip(zipped, xIn, yIn);

        ys.resize(yIn.size());
        lowess(xIn, yIn, this->GetSpanValue(), this->GetDeltaValue(), 2, ys);

        x1 = arma::rowvec(xIn);
        y1 = arma::rowvec(ys);
        // approssimo su x_out
        yOut.row(k) = util::approx(x1, y1, outputGrid);
    }

    // compute dissimilarity whit others
    arma::rowvec dso(nObs);

    for (unsigned int i = 0;i < nObs;++i)
        dso(i) = distanceObject->GetDistance(outputGrid, inputGrid.row(i), yOut, util::GetObservation(inputValues, i));

    outputCenter.Grid = outputGrid;
    outputCenter.Values = yOut;
    outputCenter.Distances = dso;

    return outputCenter;
}

CenterObject UnitQuaternionMean::GetCenter(const arma::mat& inputGrid,
                                           const arma::cube& inputValues,
                                           std::shared_ptr<Dissimilarity>& distanceObject,
                                           const arma::rowvec& outputGrid)
{
    CenterObject outputCenter;

    unsigned int nOut = outputGrid.size();
    unsigned int nDim = inputValues.n_slices;
    unsigned int nObs = inputValues.n_rows;
    unsigned int nPts = inputValues.n_cols;

    if (inputGrid.n_cols != nPts)
        Rcpp::stop("The number of columns in x should match the second dimension of y.");

    if (inputGrid.n_rows != nObs)
        Rcpp::stop("The number of rows in x should match the first dimension of y.");

    arma::mat yOut(nDim, nOut);

    // compute dissimilarity whit others
    arma::rowvec dso(nObs);
    for (unsigned int i = 0;i < nObs;++i)
        dso(i) = distanceObject->GetDistance(outputGrid, inputGrid.row(i), yOut, util::GetObservation(inputValues, i));

    outputCenter.Grid = outputGrid;
    outputCenter.Values = yOut;
    outputCenter.Distances = dso;

    return outputCenter;
}

CenterObject PseudoMedoid::GetCenter(const arma::mat& inputGrid,
                                     const arma::cube& inputValues,
                                     std::shared_ptr<Dissimilarity>& distanceObject,
                                     const arma::rowvec& outputGrid)
{
    CenterObject outputCenter;

    unsigned int nOut = outputGrid.size();
    unsigned int nDim = inputValues.n_slices;
    unsigned int nObs = inputValues.n_rows;
    unsigned int nPts = inputValues.n_cols;

    if (inputGrid.n_cols != nPts)
        Rcpp::stop("The number of columns in x should match the second dimension of y.");

    if (inputGrid.n_rows != nObs)
        Rcpp::stop("The number of rows in x should match the first dimension of y.");

    double xMin = util::GetCommonLowerBound(inputGrid);
    double xMax = util::GetCommonUpperBound(inputGrid);

    if (xMin >= xMax)
        Rcpp::stop("No common grid between curves in this group");

    arma::rowvec xCommon = arma::linspace<arma::rowvec>(xMin, xMax, nOut);

    arma::mat yOut(nDim, nOut);
    arma::mat D(nObs, nObs, arma::fill::zeros);
    arma::mat inputSlice;
    arma::colvec distVec;
    arma::rowvec workVector1, workVector2;

    for (unsigned int k = 0;k < nDim;++k)
    {
        inputSlice = inputValues.slice(k);

        for (unsigned int i = 0;i < nObs;++i)
        {
            workVector1 = util::approx(inputGrid.row(i), inputSlice.row(i), xCommon);

            for (unsigned int j = i + 1;j < nObs;++j)
            {
                workVector2 = util::approx(inputGrid.row(j), inputSlice.row(j), xCommon);
                D(i, j) = distanceObject->GetDistance(xCommon, xCommon, workVector1, workVector2);
                D(j, i) = D(i, j);
            }
        }

        distVec = arma::sum(D, 1);
        unsigned int medoidIndex =  arma::index_min(distVec);
        yOut.row(k) = util::approx(inputGrid.row(medoidIndex), inputSlice.row(medoidIndex), outputGrid);
    }

    arma::rowvec dso(nObs);

    for (unsigned int i = 0;i < nObs;++i)
        dso(i) = distanceObject->GetDistance(outputGrid, inputGrid.row(i), yOut, util::GetObservation(inputValues, i));

    outputCenter.Grid = outputGrid;
    outputCenter.Values = yOut;
    outputCenter.Distances = dso;

    return outputCenter;
}

CenterObject Median::GetCenter(const arma::mat& inputGrid,
                               const arma::cube& inputValues,
                               std::shared_ptr<Dissimilarity>& distanceObject,
                               const arma::rowvec& outputGrid)
{
    CenterObject outputCenter;

    unsigned int nOut = outputGrid.size();
    unsigned int nDim = inputValues.n_slices;
    unsigned int nObs = inputValues.n_rows;
    unsigned int nPts = inputValues.n_cols;

    if (inputGrid.n_cols != nPts)
        Rcpp::stop("The number of columns in x should match the second dimension of y.");

    if (inputGrid.n_rows != nObs)
        Rcpp::stop("The number of rows in x should match the first dimension of y.");

    double xMin = util::GetCommonLowerBound(inputGrid);
    double xMax = util::GetCommonUpperBound(inputGrid);

    if (xMin >= xMax)
        Rcpp::stop("No common grid between curves in this group");

    arma::rowvec xCommon = arma::linspace<arma::rowvec>(xMin, xMax, nOut);
    arma::mat workMatrix(nObs, nOut), yOut(nDim, nOut);

    for (unsigned int k = 0;k < nDim;++k)
    {
        // Put everyone in common grid
        for (unsigned int i = 0;i < nObs;++i)
            workMatrix.row(i) = util::approx(inputGrid.row(i), util::GetObservation(inputValues, i), xCommon);

        // Compute median
        yOut.row(k) = arma::median(workMatrix, 0);
    }

    // compute dissimilarity whit others
    arma::rowvec dso(nObs);
    arma::mat workMatrix2(nDim, nOut);
    for (unsigned int i = 0;i < nObs;++i)
    {
        for (unsigned int k = 0;k < nDim;++k)
            workMatrix2.row(k) = util::approx(inputGrid.row(i), inputValues(arma::span(i), arma::span::all, arma::span(k)), xCommon);

        dso[i] = distanceObject->GetDistance(xCommon, xCommon, yOut, workMatrix2);
    }

    outputCenter.Grid = xCommon;
    outputCenter.Values = yOut;
    outputCenter.Distances = dso;

    return outputCenter;
}

CenterObject Medoid::GetCenterParallel(const arma::mat& inputGrid,
                                       const arma::cube& inputValues,
                                       std::shared_ptr<Dissimilarity>& distanceObject,
                                       const arma::rowvec& outputGrid,
                                       unsigned int nbThreads)
{
    CenterObject outputCenter;

    unsigned int nOut = outputGrid.size();
    unsigned int nDim = inputValues.n_slices;
    unsigned int nObs = inputValues.n_rows;
    unsigned int nPts = inputValues.n_cols;

    if (inputGrid.n_cols != nPts)
        Rcpp::stop("The number of columns in x should match the second dimension of y.");

    if (inputGrid.n_rows != nObs)
        Rcpp::stop("The number of rows in x should match the first dimension of y.");

    double xMin = util::GetCommonLowerBound(inputGrid);
    double xMax = util::GetCommonUpperBound(inputGrid);

    if (xMin >= xMax)
        Rcpp::stop("No common grid between curves in this group");

    arma::rowvec xCommon = arma::linspace<arma::rowvec>(xMin, xMax, nOut);
    arma::mat yOut, workMatrix1, workMatrix2;
    arma::field<arma::rowvec> fD(nObs);

    for (unsigned int i = 0;i < nObs;++i)
        fD(i).zeros(nObs);

#ifdef _OPENMP
#pragma omp parallel for num_threads(nbThreads)
#endif

    for (unsigned int k = 1;k <= nObs * (nObs - 1) / 2;++k)
    {
        unsigned int i = std::floor((1 + std::sqrt(8 * (double)k - 7)) / 2);
        unsigned int j = k - (i - 1) * i / 2 - 1;

        workMatrix1 = util::GetObservation(inputValues, i);
        workMatrix2 = util::GetObservation(inputValues, j);

        workMatrix1 = util::approx(inputGrid.row(i), workMatrix1, xCommon);
        workMatrix2 = util::approx(inputGrid.row(j), workMatrix2, xCommon);

        fD(i)(j) = distanceObject->GetDistance(xCommon, xCommon, workMatrix1, workMatrix2);
        fD(j)(i) = fD(i)(j);
    }

    arma::colvec distVec(nObs);

    for (unsigned int i = 0;i < nObs;++i)
        distVec(i) = arma::sum(fD(i));

    unsigned int medoidIndex = arma::index_min(distVec);

    yOut = util::GetObservation(inputValues, medoidIndex);

    outputCenter.Grid = outputGrid;
    outputCenter.Values = util::approx(inputGrid.row(medoidIndex), yOut, outputGrid);
    outputCenter.Distances = fD(medoidIndex);

    return outputCenter;
}

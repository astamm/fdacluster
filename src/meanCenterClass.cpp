#include "meanCenterClass.h"

#include <squad.h>

CenterType MeanCenterMethod::GetCenter(const arma::mat& inputGrid,
                                       const arma::cube& inputValues,
                                       const std::shared_ptr<Dissimilarity>& dissimilarityPointer)
{
    CenterType outputCenter;

    unsigned int numberOfObservations = inputValues.n_rows;
    unsigned int numberOfDimensions = inputValues.n_cols;
    unsigned int numberOfPoints = inputValues.n_slices;

    if (inputGrid.n_rows != numberOfObservations)
        Rcpp::stop("The number of rows in x should match the first dimension of y.");

    if (inputGrid.n_cols != numberOfPoints)
        Rcpp::stop("The number of columns in x should match the third dimension of y.");

    // Find common grid
    double gridLowerBound = arma::max(arma::min(inputGrid, 1));
    double gridUpperBound = arma::min(arma::max(inputGrid, 1));
    arma::rowvec outputGrid = arma::linspace<arma::rowvec>(gridLowerBound, gridUpperBound, numberOfPoints);

    arma::mat meanValue(numberOfDimensions, numberOfPoints, arma::fill::zeros);
    arma::mat workMatrix;
    arma::cube yIn(numberOfObservations, numberOfDimensions, numberOfPoints);

    if (this->GetSpace() == Euclidean)
    {
        // First interpolate to common grid
        arma::cube yIn(numberOfObservations, numberOfDimensions, numberOfPoints);
        arma::vec workVector1, workVector2, tmpGrid;

        for (unsigned int i = 0;i < numberOfObservations;++i)
        {
            tmpGrid = inputGrid.row(i);

            for (unsigned int j = 0;j < numberOfDimensions;++j)
            {
                workVector1 = inputValues(arma::span(i), arma::span(j), arma::span::all);
                arma::interp1(tmpGrid, workVector1, outputGrid, workVector2, "*linear");
                yIn(arma::span(i), arma::span(j), arma::span::all) = workVector2;
            }
        }

        // then, compute mean pointwise
        for (unsigned int i = 0;i < numberOfPoints;++i)
            meanValue.col(i) = arma::mean(yIn.slice(i), 0);
    }
    else if (this->GetSpace() == UnitQuaternion)
    {
        // First interpolate to common grid
        Rcpp::NumericVector tmpVec;
        Rcpp::NumericMatrix tmpMat;

        for (unsigned int i = 0;i < numberOfObservations;++i)
        {
            tmpVec = Rcpp::wrap(inputGrid.row(i));
            workMatrix = inputValues(arma::span(i), arma::span::all, arma::span::all);
            tmpMat = Rcpp::wrap(workMatrix);
            tmpMat = squad::RegularizeGrid(tmpVec, tmpMat, gridLowerBound, gridUpperBound, numberOfPoints);
            yIn.row(i) = Rcpp::as<arma::mat>(tmpMat);
        }

        // then, compute mean pointwise
        for (unsigned int i = 0;i < numberOfPoints;++i)
        {
            tmpMat = squad::GetGeodesicMean(Rcpp::wrap(yIn.slice(i)));
            meanValue.col(i) = Rcpp::as<arma::vec>(tmpMat);
        }
    }
    else
        Rcpp::Rcout << "Operations for the requested space are not yet implemented." << std::endl;

    // compute dissimilarity between observations and center
    arma::rowvec distancesToCenter(numberOfObservations);
    for (unsigned int i = 0;i < numberOfObservations;++i)
    {
        workMatrix = inputValues(arma::span(i), arma::span::all, arma::span::all);
        distancesToCenter(i) = dissimilarityPointer->GetDistance(
            outputGrid,
            inputGrid.row(i),
            meanValue,
            workMatrix
        );
    }

    outputCenter.centerGrid = outputGrid;
    outputCenter.centerValues = meanValue;
    outputCenter.distancesToCenter = distancesToCenter;

    return outputCenter;
}

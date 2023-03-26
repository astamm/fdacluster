#include "meanCenterClass.h"

CenterType MeanCenterMethod::GetCenter(const arma::mat& inputGrid,
                                       const arma::cube& inputValues,
                                       const std::shared_ptr<BaseDissimilarityFunction>& dissimilarityPointer)
{
    CenterType outputCenter;

    unsigned int numberOfObservations = inputValues.n_rows;
    unsigned int numberOfDimensions = inputValues.n_cols;
    unsigned int numberOfPoints = inputValues.n_slices;

    if (inputGrid.n_rows != numberOfObservations)
        Rcpp::stop("The number of rows in x should match the first dimension of y.");

    if (inputGrid.n_cols != numberOfPoints)
        Rcpp::stop("The number of columns in x should match the third dimension of y.");

    // Find largest grid
    double gridLowerBound = inputGrid.min();
    double gridUpperBound = inputGrid.max();
    arma::rowvec outGrid = arma::linspace<arma::rowvec>(gridLowerBound, gridUpperBound, numberOfPoints);

    arma::uvec finiteIndices;
    arma::mat meanValue(numberOfDimensions, numberOfPoints, arma::fill::zeros);
    arma::mat workMatrix;
    arma::cube yIn(numberOfObservations, numberOfDimensions, numberOfPoints);

    // First interpolate to common grid
    arma::rowvec inGrid;
    arma::rowvec inValue;
    arma::rowvec outValue;

    for (unsigned int i = 0;i < numberOfObservations;++i)
    {
      inGrid = inputGrid.row(i);

      for (unsigned int j = 0;j < numberOfDimensions;++j)
      {
        inValue = inputValues.tube(i, j);
        arma::interp1(inGrid, inValue, outGrid, outValue, "*linear");
        yIn.tube(i, j) = outValue;
      }
    }

    // Next, compute point-wise mean
    for (unsigned int i = 0;i < numberOfPoints;++i)
    {
      finiteIndices = arma::find_finite(yIn.slice(i).col(0));
      workMatrix = yIn.slice(i).rows(finiteIndices);
      meanValue.col(i) = arma::mean(workMatrix, 0).as_col();
    }

    // Finally, compute dissimilarity between observations and center
    arma::rowvec distancesToCenter(numberOfObservations);
    for (unsigned int i = 0;i < numberOfObservations;++i)
    {
        workMatrix = inputValues.row(i);
        distancesToCenter(i) = dissimilarityPointer->GetDistance(
            outGrid,
            inputGrid.row(i),
            meanValue,
            workMatrix
        );
    }

    outputCenter.centerGrid = outGrid;
    outputCenter.centerValues = meanValue;
    outputCenter.distancesToCenter = distancesToCenter;

    return outputCenter;
}

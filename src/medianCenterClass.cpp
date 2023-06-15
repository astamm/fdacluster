#include "medianCenterClass.h"

CenterType MedianCenterMethod::GetCenter(const arma::mat& inputGrid,
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

  // Find intersection grid
  double gridLowerBound = inputGrid.col(0).max();
  double gridUpperBound = inputGrid.col(numberOfPoints - 1).min();
  arma::rowvec outGrid = arma::linspace<arma::rowvec>(gridLowerBound, gridUpperBound, numberOfPoints);

  arma::uvec finiteIndices;
  arma::mat medianValue(numberOfDimensions, numberOfPoints, arma::fill::zeros);
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
    medianValue.col(i) = arma::median(workMatrix, 0).as_col();
  }

  // Finally, compute dissimilarity between observations and center
  arma::rowvec distancesToCenter(numberOfObservations);
  for (unsigned int i = 0;i < numberOfObservations;++i)
  {
    workMatrix = inputValues.row(i);
    distancesToCenter(i) = dissimilarityPointer->GetDistance(
      outGrid,
      inputGrid.row(i),
      medianValue,
      workMatrix
    );
  }

  outputCenter.centerGrid = outGrid;
  outputCenter.centerValues = medianValue;
  outputCenter.distancesToCenter = distancesToCenter;

  return outputCenter;
}

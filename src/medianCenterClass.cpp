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
  // double gridLowerBound = inputGrid.col(0).max();
  // double gridUpperBound = inputGrid.col(numberOfPoints - 1).min();
  // arma::rowvec outGrid = arma::linspace<arma::rowvec>(gridLowerBound, gridUpperBound, numberOfPoints);

  // Define union grid
  double unionLowerBound = inputGrid.min();
  double unionUpperBound = inputGrid.max();
  arma::rowvec outGrid = arma::linspace<arma::rowvec>(unionLowerBound, unionUpperBound, numberOfPoints);

  arma::uvec finiteIndices, nonFiniteIndices;
  arma::mat centroidValue(numberOfDimensions, numberOfPoints, arma::fill::zeros);
  arma::mat workMatrix(numberOfObservations, numberOfDimensions);
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
  arma::uvec numberOfMissingPoints(numberOfPoints);
  for (unsigned int i = 0;i < numberOfPoints;++i)
  {
    finiteIndices = arma::find_finite(yIn.slice(i).col(0));
    numberOfMissingPoints(i) = numberOfObservations - finiteIndices.size();
    workMatrix = yIn.slice(i).rows(finiteIndices);
    centroidValue.col(i) = arma::median(workMatrix, 0).as_col();
  }

  // Next, compute point-wise median
  int halfPoint = (int)(numberOfPoints / 2);
  arma::rowvec localSlopeValues, localInterceptValues, localPredictedValues;
  arma::rowvec shiftValues;

  for (int i = halfPoint;i < numberOfPoints;++i)
  {
    nonFiniteIndices = arma::find_nonfinite(yIn.slice(i).col(0));
    for (unsigned int j = 0;j < nonFiniteIndices.size();++j)
    {
      unsigned int rowIndex = nonFiniteIndices(j);
      localSlopeValues = (yIn.slice(i - 1).row(rowIndex) - yIn.slice(i - 2).row(rowIndex)) / (outGrid(i - 1) - outGrid(i - 2));
      localInterceptValues = yIn.slice(i - 2).row(rowIndex) - localSlopeValues * outGrid(i - 2);
      localPredictedValues = localInterceptValues + localSlopeValues * outGrid(i);
      shiftValues = localPredictedValues - centroidValue.col(i).as_row();
      for (int k = i;k < numberOfPoints;++k)
      {
        if (k > i && numberOfMissingPoints(k) != numberOfMissingPoints(k - 1))
          break;
        yIn.slice(k).row(rowIndex) = centroidValue.col(k).as_row() + shiftValues;
      }
    }

    centroidValue.col(i) = arma::median(yIn.slice(i), 0).as_col();
  }

  for (int i = halfPoint - 1;i >= 0;--i)
  {
    nonFiniteIndices = arma::find_nonfinite(yIn.slice(i).col(0));
    for (unsigned int j = 0;j < nonFiniteIndices.size();++j)
    {
      unsigned int rowIndex = nonFiniteIndices(j);
      localSlopeValues = (yIn.slice(i + 1).row(rowIndex) - yIn.slice(i + 2).row(rowIndex)) / (outGrid(i + 1) - outGrid(i + 2));
      localInterceptValues = yIn.slice(i + 2).row(rowIndex) - localSlopeValues * outGrid(i + 2);
      localPredictedValues = localInterceptValues + localSlopeValues * outGrid(i);
      shiftValues = localPredictedValues - centroidValue.col(i).as_row();
      for (int k = i;k >= 0;--k)
      {
        if (k < i && numberOfMissingPoints(k) != numberOfMissingPoints(k + 1))
          break;
        yIn.slice(k).row(rowIndex) = centroidValue.col(k).as_row() + shiftValues;
      }
    }

    centroidValue.col(i) = arma::median(yIn.slice(i), 0).as_col();
  }

  // Finally, compute dissimilarity between observations and center
  arma::rowvec distancesToCenter(numberOfObservations);
  for (unsigned int i = 0;i < numberOfObservations;++i)
  {
    workMatrix = inputValues.row(i);
    distancesToCenter(i) = dissimilarityPointer->GetDistance(
      outGrid,
      inputGrid.row(i),
      centroidValue,
      workMatrix
    );
  }

  outputCenter.centerGrid = outGrid;
  outputCenter.centerValues = centroidValue;
  outputCenter.distancesToCenter = distancesToCenter;

  return outputCenter;
}

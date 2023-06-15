#include "polyCenterClass.h"

CenterType PolyCenterMethod::GetCenter(const arma::mat& inputGrid,
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

  arma::mat meanValue(numberOfDimensions, numberOfPoints, arma::fill::zeros);
  arma::mat workMatrix;
  std::vector<double> xIn, yIn;
  arma::rowvec inGrid, inValue, outValue(outGrid.size());
  arma::vec polyCoefficients(this->GetPolynomialDegree() + 1);

  for (unsigned int l = 0;l < numberOfDimensions;++l)
  {
    xIn.clear();
    yIn.clear();

    for (unsigned int i = 0;i < numberOfObservations;++i)
    {
      for (unsigned int j = 0;j < numberOfPoints;++j)
      {
        if (arma::is_finite(inputGrid(i, j)) && arma::is_finite(inputValues(i, l, j)))
        {
          xIn.push_back(inputGrid(i, j));
          yIn.push_back(inputValues(i, l, j));
        }
      }
    }

    inGrid = arma::conv_to<arma::rowvec>::from(xIn);
    inValue = arma::conv_to<arma::rowvec>::from(yIn);
    arma::polyfit(polyCoefficients, inGrid, inValue, this->GetPolynomialDegree());
    outValue = arma::polyval(polyCoefficients, outGrid);
    meanValue.row(l) = outValue;
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

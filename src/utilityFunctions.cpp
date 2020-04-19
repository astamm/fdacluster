#include "utilityFunctions.h"

#include <squad.h>

//
// tableCpp
//
std::map<unsigned int, unsigned int> tableCpp(const arma::urowvec &inputLabels)
{
    std::map<unsigned int, unsigned int> outputCounts;
    unsigned int numberOfObservations = inputLabels.size();

    for (unsigned int i = 0;i < numberOfObservations;++i)
        ++outputCounts[inputLabels[i]];

    return outputCounts;
}

//
//  which_out
//
arma::urowvec which_out(const arma::rowvec &inputValues, const double lowerBound, const double upperBound)
{
    unsigned int nbValues = inputValues.size();
    arma::urowvec outputIndices(nbValues);
    unsigned int pos = 0;

    for (unsigned int i = 0;i < nbValues;++i)
    {
        if (inputValues(i) < lowerBound || inputValues(i) > upperBound)
        {
            outputIndices(pos) = i;
            ++pos;
        }
    }

    outputIndices.resize(pos);

    return outputIndices;
}

//
//  quantile
//
double quantile(const arma::rowvec &inputValues, const double quantileOrder)
{
    if (quantileOrder <= 0 || quantileOrder >= 1)
        Rcpp::stop("The quantile order should be between 0 and 1 excluded.");

    std::vector<double> workVector = arma::conv_to<std::vector<double> >::from(inputValues);
    unsigned int q = std::round(inputValues.size() * quantileOrder);

    std::nth_element(workVector.begin(), workVector.begin() + q, workVector.end());

    return *(workVector.begin() + q);
}

//
// GetCommonLowerBound
//
double GetCommonLowerBound(const arma::mat& inputGrids)
{
    // Assumption: Matrix is in format NOBS x NPTS
    unsigned int nObs = inputGrids.n_rows;
    unsigned int nPts = inputGrids.n_cols;

    double resVal = 0.0;
    unsigned int posi = 0;

    for (unsigned int i = 0;i < nObs;++i)
    {
        double tmpVal = 0.0;
        unsigned int posj = 0;

        for (unsigned int j = 0;j < nPts;++j)
        {
            if (!arma::is_finite(inputGrids(i, j)))
                continue;

            if (inputGrids(i, j) < tmpVal || posj == 0)
                tmpVal = inputGrids(i, j);

            ++posj;
        }

        if (tmpVal > resVal || posi == 0)
            resVal = tmpVal;

        ++posi;
    }

    return resVal;
}

//
// GetCommonUpperBound
//
double GetCommonUpperBound(const arma::mat& inputGrids)
{
    // Assumption: Matrix is in format NOBS x NPTS
    unsigned int nObs = inputGrids.n_rows;
    unsigned int nPts = inputGrids.n_cols;

    double resVal = 0.0;
    unsigned int posi = 0;

    for (unsigned int i = 0;i < nObs;++i)
    {
        double tmpVal = 0.0;
        unsigned int posj = 0;

        for (unsigned int j = 0;j < nPts;++j)
        {
            if (!arma::is_finite(inputGrids(i, j)))
                continue;

            if (inputGrids(i, j) > tmpVal || posj == 0)
                tmpVal = inputGrids(i, j);

            ++posj;
        }

        if (tmpVal < resVal || posi == 0)
            resVal = tmpVal;

        ++posi;
    }

    return resVal;
}

//
// observations
//
arma::cube GetObservations(const arma::cube& inputData, arma::uvec& observationIndices)
{
    arma::cube outputCube(observationIndices.size(), inputData.n_cols, inputData.n_slices);

    for (unsigned int i = 0;i < observationIndices.size();++i)
        outputCube.row(i) = inputData.row(observationIndices(i));

    return outputCube;
}

//
// approx
//
Rcpp::List approx(const arma::rowvec& inputGrid,
                  const arma::mat& inputValues,
                  const std::string interpolationMethod)
{
    // inputGrid is assumed to be of size NDIM x NPTS
    // inputValues is assumed to be of size NPTS
    // outputX is assumed to be of size NOUT
    // outputY will be of size NDIM x NPTS
    unsigned int nPts = inputValues.n_cols;
    unsigned int nDim = inputValues.n_rows;

    if (inputGrid.size() != nPts)
        Rcpp::stop("The length of input arc length and the number of rows in the input matrix should be equal.");

    arma::rowvec inputXCopy(inputGrid);
    arma::mat inputYCopy(inputValues);

    for (int j = nPts - 1;j >= 0;--j)
    {
        if (arma::is_finite(inputGrid(j)) && arma::is_finite(inputValues.col(j)))
            continue;

        inputXCopy.shed_col(j);
        inputYCopy.shed_col(j);
    }

    nPts = inputXCopy.size();

    if (nPts <= 1)
        Rcpp::stop("All entries are NA except maybe one. Interpolation is not possible. Aborting...");

    double gridMin = inputXCopy.min();
    double gridMax = inputXCopy.max();
    arma::rowvec outputGrid = arma::linspace<arma::rowvec>(gridMin, gridMax, nPts);

    arma::mat outputValues(nPts, nDim);
    outputValues.fill(arma::datum::nan);

    if (outputGrid(nPts - 1) - outputGrid(0) <= 0)
    {
        Rcpp::Rcout << "Output arc length is not in increasing order: " << outputGrid(nPts - 1) - outputGrid(0) << std::endl;
        return Rcpp::List::create(
            Rcpp::Named("grid") = outputGrid,
            Rcpp::Named("values") = outputValues
        );
    }

    if (interpolationMethod == "Linear")
    {
        unsigned int i = 0;

        while (i < nPts && outputGrid(i) <= inputXCopy(0))
        {
            if (std::abs(outputGrid(i) - inputXCopy(0)) < arma::datum::eps)
                outputValues.col(i) = inputYCopy.col(0);
            ++i;
        }

        if (i == nPts)
            return Rcpp::List::create(
                Rcpp::Named("grid") = outputGrid,
                Rcpp::Named("values") = outputValues
            );

        arma::vec oldY, newY;

        for (unsigned int j = 1;j < nPts;++j)
        {
            double oldX = inputXCopy(j - 1);
            oldY = inputYCopy.col(j - 1);
            double newX = inputXCopy(j);
            newY = inputYCopy.col(j);

            while (i < nPts && outputGrid(i) <= newX)
            {
                outputValues.col(i) = oldY + (newY - oldY) * (outputGrid(i) - oldX) / (newX - oldX);
                ++i;
            }
        }
    }
    else if (interpolationMethod == "UnitQuaternion")
    {
        Rcpp::NumericVector workingGrid = Rcpp::wrap(inputXCopy);
        Rcpp::NumericMatrix workingValues = Rcpp::wrap(inputYCopy);
        workingValues = squad::RegularizeGrid(workingGrid, workingValues, gridMin, gridMax, nPts);
        inputYCopy = Rcpp::as<arma::mat>(workingValues);
    }
    else
        Rcpp::stop("Unsupported interpolation method. Currently supported ones are linear or unit quaternion.");

    return Rcpp::List::create(
        Rcpp::Named("grid") = outputGrid,
        Rcpp::Named("values") = outputValues
    );
}

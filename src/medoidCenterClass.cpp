#include "medoidCenterClass.h"

#ifdef _OPENMP
#include <omp.h>
#endif

CenterType MedoidCenterMethod::GetCenter(const arma::mat& inputGrid,
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

    arma::mat distanceMatrix(numberOfObservations, numberOfObservations, arma::fill::zeros);
    // These auxiliary matrices will be of size nDim x nPts
    arma::mat workMatrix1, workMatrix2;

    for (unsigned int i = 0;i < numberOfObservations;++i)
    {
        workMatrix1 = inputValues(arma::span(i), arma::span::all, arma::span::all);

        for (unsigned int j = i + 1;j < numberOfObservations;++j)
        {
            workMatrix2 = inputValues(arma::span(j), arma::span::all, arma::span::all);

            double workDistance = dissimilarityPointer->GetDistance(
                inputGrid.row(i),
                inputGrid.row(j),
                workMatrix1,
                workMatrix2
            );

            distanceMatrix(i, j) = workDistance;
            distanceMatrix(j, i) = workDistance;
        }
    }

    arma::colvec distanceVector = arma::sum(distanceMatrix, 1);
    unsigned int medoidIndex = arma::index_min(distanceVector);

    outputCenter.centerGrid = inputGrid.row(medoidIndex);
    outputCenter.centerValues = inputValues(arma::span(medoidIndex), arma::span::all, arma::span::all);
    outputCenter.distancesToCenter = distanceVector.row(medoidIndex);

    return outputCenter;
}

CenterType MedoidCenterMethod::GetCenterParallel(const arma::mat& inputGrid,
                                                 const arma::cube& inputValues,
                                                 const std::shared_ptr<Dissimilarity>& dissimilarityPointer,
                                                 unsigned int nbThreads)
{
    CenterType outputCenter;

    unsigned int numberOfObservations = inputValues.n_rows;
    unsigned int numberOfDimensions = inputValues.n_cols;
    unsigned int numberOfPoints = inputValues.n_slices;

    if (inputGrid.n_rows != numberOfObservations)
        Rcpp::stop("The number of rows in x should match the first dimension of y.");

    if (inputGrid.n_cols != numberOfPoints)
        Rcpp::stop("The number of columns in x should match the third dimension of y.");

    arma::mat distanceMatrix(numberOfObservations, numberOfObservations, arma::fill::zeros);
    // These auxiliary matrices will be of size nDim x nPts
    arma::mat workMatrix1, workMatrix2;

    arma::field<arma::rowvec> distanceField(numberOfObservations);
    for (unsigned int i = 0;i < numberOfObservations;++i)
        distanceField(i).zeros(numberOfObservations);

#ifdef _OPENMP
#pragma omp parallel for num_threads(nbThreads)
#endif

    for (unsigned int k = 1;k <= numberOfObservations * (numberOfObservations - 1) / 2;++k)
    {
        unsigned int i = std::floor((1 + std::sqrt(8 * (double)k - 7)) / 2);
        unsigned int j = k - (i - 1) * i / 2 - 1;

        workMatrix1 = inputValues(arma::span(i), arma::span::all, arma::span::all);
        workMatrix2 = inputValues(arma::span(j), arma::span::all, arma::span::all);

        double workDistance = dissimilarityPointer->GetDistance(
            inputGrid.row(i),
            inputGrid.row(j),
            workMatrix1,
            workMatrix2
        );

        distanceMatrix(i, j) = workDistance;
        distanceMatrix(j, i) = workDistance;
    }

    arma::colvec distanceVector = arma::sum(distanceMatrix, 1);
    unsigned int medoidIndex = arma::index_min(distanceVector);

    outputCenter.centerGrid = inputGrid.row(medoidIndex);
    outputCenter.centerValues = inputValues(arma::span(medoidIndex), arma::span::all, arma::span::all);
    outputCenter.distancesToCenter = distanceVector.row(medoidIndex);

    return outputCenter;
}

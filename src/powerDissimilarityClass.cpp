#include "powerDissimilarityClass.h"

#include <squad.h>

double PowerDissimilarityFunction::GetDistance(const arma::rowvec& grid1,
                                               const arma::rowvec& grid2,
                                               const arma::mat& values1,
                                               const arma::mat& values2)
{
    FunctionPairType pair = this->GetComparableFunctions(grid1, grid2, values1, values2);

    if (pair.Grid.is_empty())
        return DBL_MAX;

    unsigned int nDim = pair.Values1.n_rows;
    unsigned int nPts = pair.Grid.size();

    if (nPts <= 1.0)
        return DBL_MAX;

    double squaredDistanceValue = 0.0;
    double squaredNorm1Value = 0.0;
    double squaredNorm2Value = 0.0;

    arma::rowvec weightVector(nPts - 1, arma::fill::zeros);
    for (unsigned int i = 0;i < nPts - 1;++i)
    {
        double denomValue = (m_DomainUpperBound - m_DomainLowerBound) * (i + 1.0) / (nPts - 1.0);

        if (denomValue < std::sqrt(std::numeric_limits<double>::epsilon()))
            continue;

        weightVector[i] = (m_DomainUpperBound - m_DomainLowerBound) / denomValue;
    }

    if (this->GetSpace() == Euclidean)
    {
        arma::rowvec workVector;

        for (unsigned int k = 0;k < nDim;++k)
        {
            workVector = pair.Values1.row(k).cols(1, nPts - 1) - pair.Values2.row(k).cols(1, nPts - 1);
            squaredDistanceValue += arma::dot(workVector, weightVector % workVector);
            workVector = pair.Values1.row(k).cols(1, nPts - 1);
            squaredNorm1Value += arma::dot(workVector, weightVector % workVector);
            workVector = pair.Values2.row(k).cols(1, nPts - 1);
            squaredNorm2Value += arma::dot(workVector, weightVector % workVector);
        }
    }
    else if (this->GetSpace() == UnitQuaternion)
    {
        arma::mat neutralMatrix(nDim, nPts, arma::fill::zeros);
        neutralMatrix.row(0).ones();

        for (unsigned int j = 0;j < nPts - 1;++j)
        {
            double workValue = squad::GeodesicQuaternionDistance(
                Rcpp::wrap(pair.Values1),
                Rcpp::wrap(pair.Values2),
                j + 1, j + 1
            );

            squaredDistanceValue += weightVector[j] * workValue * workValue;

            workValue = squad::GeodesicQuaternionDistance(
                Rcpp::wrap(pair.Values1),
                Rcpp::wrap(neutralMatrix),
                j + 1, j + 1
            );

            squaredNorm1Value += weightVector[j] * workValue * workValue;

            workValue = squad::GeodesicQuaternionDistance(
                Rcpp::wrap(pair.Values2),
                Rcpp::wrap(neutralMatrix),
                j + 1, j + 1
            );

            squaredNorm2Value += weightVector[j] * workValue * workValue;
        }
    }
    else
        Rcpp::Rcout << "Distance operations for the requested space are not yet implemented." << std::endl;

    return std::sqrt(squaredDistanceValue) / (std::sqrt(squaredNorm1Value) + std::sqrt(squaredNorm2Value));
}

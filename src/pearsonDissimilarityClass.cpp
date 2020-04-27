#include "pearsonDissimilarityClass.h"

#include <squad.h>

double PearsonDissimilarityFunction::GetDistance(const arma::rowvec& grid1,
                                                 const arma::rowvec& grid2,
                                                 const arma::mat& values1,
                                                 const arma::mat& values2)
{
    FunctionPairType pair = this->GetComparableFunctions(grid1, grid2, values1, values2);

    if (pair.Grid.is_empty())
        return DBL_MAX;

    unsigned int nDim = pair.Values1.n_rows;
    unsigned int nPts = pair.Values1.n_cols;

    double squaredDistanceValue = 0.0;
    double squaredNorm1Value = 0.0;
    double squaredNorm2Value = 0.0;

    if (this->GetSpace() == Euclidean)
    {
        arma::rowvec workVector;

        for (unsigned int k = 0;k < nDim;++k)
        {
            workVector = pair.Values1.row(k).cols(1, nPts - 1);
            squaredNorm1Value += arma::dot(workVector, workVector);
            workVector = pair.Values2.row(k).cols(1, nPts - 1);
            squaredNorm2Value += arma::dot(workVector, workVector);
        }

        squaredNorm1Value /= ((double)nPts - 1.0);
        squaredNorm2Value /= ((double)nPts - 1.0);

        for (unsigned int k = 0;k < nDim;++k)
        {
            workVector  = pair.Values1.row(k).cols(1, nPts - 1) / std::sqrt(squaredNorm1Value);
            workVector -= pair.Values2.row(k).cols(1, nPts - 1) / std::sqrt(squaredNorm2Value);
            squaredDistanceValue += arma::dot(workVector, workVector);
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
                Rcpp::wrap(neutralMatrix),
                j + 1, j + 1
            );

            squaredNorm1Value += workValue * workValue;

            workValue = squad::GeodesicQuaternionDistance(
                Rcpp::wrap(pair.Values2),
                Rcpp::wrap(neutralMatrix),
                j + 1, j + 1
            );

            squaredNorm2Value += workValue * workValue;
        }

        squaredNorm1Value /= ((double)nPts - 1.0);
        squaredNorm2Value /= ((double)nPts - 1.0);

        for (unsigned int j = 0;j < nPts - 1;++j)
        {
            double workValue = squad::GeodesicQuaternionDistance(
                Rcpp::wrap(pair.Values1 / std::sqrt(squaredNorm1Value)),
                Rcpp::wrap(pair.Values2 / std::sqrt(squaredNorm2Value)),
                j + 1, j + 1
            );

            squaredDistanceValue += workValue * workValue;
        }
    }
    else
        Rcpp::Rcout << "Distance operations for the requested space are not yet implemented." << std::endl;

    squaredDistanceValue /= ((double)nPts - 1.0);

    return std::sqrt(squaredDistanceValue) / 2.0;
}

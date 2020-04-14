#include "l2DissimilarityClass.h"

#include <squad.h>

double L2DissimilarityFunction::GetDistance(const arma::rowvec& grid1,
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

    arma::rowvec diffVector = pair.Grid.cols(1, nPts - 1) - pair.Grid.cols(0, nPts - 2);
    double rangeValue = pair.Grid(nPts - 1) - pair.Grid(0);
    double squaredDistanceValue = 0.0;

    if (this->GetSpace() == Euclidean)
    {
        arma::rowvec workVector;

        for (unsigned int k = 0;k < nDim;++k)
        {
            workVector = arma::sqrt(diffVector) % (pair.Values1.row(k).cols(1, nPts - 1) - pair.Values2.row(k).cols(1, nPts - 1));
            squaredDistanceValue += arma::dot(workVector, workVector);
        }

        squaredDistanceValue /= (rangeValue * (double)nDim);
    }
    else if (this->GetSpace() == UnitQuaternion)
    {
        for (unsigned int j = 0;j < nPts - 1;++j)
        {
            double tmpDistance = squad::GeodesicQuaternionDistance(
                Rcpp::wrap(pair.Values1),
                Rcpp::wrap(pair.Values2),
                j + 1, j + 1
            );

            squaredDistanceValue += diffVector[j] * tmpDistance * tmpDistance;
        }

        squaredDistanceValue /= rangeValue;
    }
    else
        Rcpp::Rcout << "Distance operations for the requested space are not yet implemented." << std::endl;

    return std::sqrt(squaredDistanceValue);
}

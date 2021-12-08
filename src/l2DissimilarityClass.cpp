#include "l2DissimilarityClass.h"

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

    double squaredDistanceValue = 0.0;
    double squaredNorm1Value = 0.0;
    double squaredNorm2Value = 0.0;

    arma::rowvec workVector;

    for (unsigned int k = 0;k < nDim;++k)
    {
      workVector = pair.Values1.row(k).cols(1, nPts - 1) - pair.Values2.row(k).cols(1, nPts - 1);
      squaredDistanceValue += arma::dot(workVector, workVector);
      workVector = pair.Values1.row(k).cols(1, nPts - 1);
      squaredNorm1Value += arma::dot(workVector, workVector);
      workVector = pair.Values2.row(k).cols(1, nPts - 1);
      squaredNorm2Value += arma::dot(workVector, workVector);
    }

    return std::sqrt(squaredDistanceValue) / (std::sqrt(squaredNorm1Value) + std::sqrt(squaredNorm2Value));
}

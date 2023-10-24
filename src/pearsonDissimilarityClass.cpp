#include "pearsonDissimilarityClass.h"

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

    arma::rowvec workVector;

    for (unsigned int k = 0;k < nDim;++k)
    {
      workVector = pair.Values1.row(k).cols(1, nPts - 1);
      squaredNorm1Value += arma::dot(workVector, workVector);
      workVector = pair.Values2.row(k).cols(1, nPts - 1);
      squaredNorm2Value += arma::dot(workVector, workVector);
    }

    double epsValue = std::sqrt(std::numeric_limits<double>::epsilon());
    if (squaredNorm1Value < epsValue && squaredNorm2Value < epsValue)
      return 0.0;

    squaredNorm1Value /= ((double)nPts - 1.0);
    squaredNorm2Value /= ((double)nPts - 1.0);

    for (unsigned int k = 0;k < nDim;++k)
    {
      workVector  = pair.Values1.row(k).cols(1, nPts - 1) / std::sqrt(squaredNorm1Value);
      workVector -= pair.Values2.row(k).cols(1, nPts - 1) / std::sqrt(squaredNorm2Value);
      squaredDistanceValue += arma::dot(workVector, workVector);
    }

    squaredDistanceValue /= ((double)nPts - 1.0);

    return std::sqrt(squaredDistanceValue) / 2.0;
}

#include "normalizedL2DissimilarityClass.h"

double NormalizedL2DissimilarityFunction::GetDistance(const arma::rowvec& grid1,
                                                      const arma::rowvec& grid2,
                                                      const arma::mat& values1,
                                                      const arma::mat& values2)
{
    FunctionPairType pair = this->GetComparableFunctions(grid1, grid2, values1, values2);

    if (pair.Grid.is_empty())
        return DBL_MAX;

    unsigned int nPts = pair.Grid.size();

    if (nPts <= 1.0)
        return DBL_MAX;

    double squaredNorm1Value = arma::sum(arma::trapz(pair.Grid, arma::pow(pair.Values1, 2.0), 1));
    double squaredNorm2Value = arma::sum(arma::trapz(pair.Grid, arma::pow(pair.Values2, 2.0), 1));

    double epsValue = std::sqrt(std::numeric_limits<double>::epsilon());
    if (squaredNorm1Value < epsValue && squaredNorm2Value < epsValue)
      return 0.0;

    double squaredDistanceValue = arma::sum(arma::trapz(pair.Grid, arma::pow(pair.Values1 - pair.Values2, 2.0), 1));

    return std::sqrt(squaredDistanceValue) / (std::sqrt(squaredNorm1Value) + std::sqrt(squaredNorm2Value));
}

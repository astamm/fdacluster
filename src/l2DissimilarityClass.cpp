#include "l2DissimilarityClass.h"

double L2DissimilarityFunction::GetDistance(const arma::rowvec& grid1,
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

    return std::sqrt(arma::sum(arma::trapz(pair.Grid, arma::pow(pair.Values1 - pair.Values2, 2.0), 1)));
}

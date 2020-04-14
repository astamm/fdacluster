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
    double pearsonCorrelation = 0.0;

    for (unsigned int k = 0;k < nDim;++k)
        pearsonCorrelation += arma::norm_dot(pair.Values1.row(k), pair.Values2.row(k));

    return -pearsonCorrelation / (double)nDim;
}

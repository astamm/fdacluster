#include "baseWarpingClass.h"

WarpingSet BaseWarpingFunction::SetInputData(const arma::rowvec &grid1,
                                             const arma::rowvec &grid2,
                                             const arma::mat &values1,
                                             const arma::mat &values2,
                                             const std::shared_ptr<Dissimilarity> &dissimilarity)
{
    WarpingSet out;

    out.inputGrid1 = grid1;
    out.inputGrid2 = grid2;
    out.inputValues1 = values1;
    out.inputValues2 = values2;
    out.dissimilarityPointer = dissimilarity;

    return out;
}

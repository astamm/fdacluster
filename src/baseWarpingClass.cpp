#include "baseWarpingClass.h"

WarpingSet BaseWarpingFunction::SetInputData(const arma::rowvec &grid1,
                                             const arma::rowvec &grid2,
                                             const arma::mat &values1,
                                             const arma::mat &values2,
                                             const std::shared_ptr<BaseDissimilarityFunction> &dissimilarity)
{
    WarpingSet out;

    out.inputGrid1 = grid1;
    out.inputGrid2 = grid2;
    out.inputValues1 = values1;
    out.inputValues2 = values2;
    out.dissimilarityPointer = dissimilarity;

    return out;
}

void BaseWarpingFunction::SetParameterBounds(const arma::mat &bounds)
{
    m_ParameterLowerBounds = bounds.row(0);
    m_ParameterUpperBounds = bounds.row(1);
}

double BaseWarpingFunction::GetDissimilarityAfterWarping(const WarpingSet &warpingSet,
                                                         const arma::rowvec &warpingParameters)
{
    return warpingSet.dissimilarityPointer->GetDistance(
            this->ApplyWarping(warpingSet.inputGrid1, warpingParameters),
            warpingSet.inputGrid2,
            warpingSet.inputValues1,
            warpingSet.inputValues2
    );
}

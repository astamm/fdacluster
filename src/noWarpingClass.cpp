#include "noWarpingClass.h"

unsigned int NoWarpingFunction::GetNumberOfParameters()
{
    return 0;
}

arma::mat NoWarpingFunction::ApplyWarping(const arma::mat &inputGrids,
                                          const arma::mat &warpingParameters)
{
    return inputGrids;
}

void NoWarpingFunction::SetParameterBounds(const arma::rowvec &warpingOptions,
                                           const arma::mat &inputGrids)
{
    m_ParameterLowerBounds.set_size(0);
    m_ParameterUpperBounds.set_size(0);
}

arma::mat NoWarpingFunction::GetFinalWarping(const arma::cube &warpingParametersContainer,
                                             const arma::urowvec &observationMemberships,
                                             const arma::urowvec &clusterIndices)
{
    arma::mat out(observationMemberships.n_cols, 0);
    return out;
}

void NoWarpingFunction::Normalize(arma::mat &warpingParameters,
                                  const arma::urowvec &clusterIndices,
                                  const arma::urowvec &observationMemberships)
{
    return;
}

double NoWarpingFunction::GetDissimilarityAfterWarping(const WarpingSet &warpingSet,
                                                       const arma::rowvec &warpingParameters)
{
    return warpingSet.dissimilarityPointer->GetDistance(
            warpingSet.inputGrid1,
            warpingSet.inputGrid2,
            warpingSet.inputValues1,
            warpingSet.inputValues2
    );
}

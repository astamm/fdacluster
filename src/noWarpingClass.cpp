#include "noWarpingClass.h"

arma::mat NoWarpingFunction::ApplyWarping(const arma::mat &x, const arma::mat &par)
{
    return x;
}

void NoWarpingFunction::SetParameterBounds(const arma::rowvec &war_opt,
                                           const arma::mat &x)
{
    m_ParameterLowerBounds.set_size(0);
    m_ParameterUpperBounds.set_size(0);
}

arma::mat NoWarpingFunction::GetFinalWarping(const arma::cube &parameters_vec,
                                             const arma::urowvec &labels,
                                             const arma::urowvec &ict)
{
    arma::mat out(0, labels.n_cols);
    return out;
}

void NoWarpingFunction::Normalize(arma::mat &par,
                                  const arma::urowvec &ict,
                                  const arma::urowvec &labels)
{
    return;
}

double NoWarpingFunction::GetDissimilarityAfterWarping(const WarpingSet &warpingSet,
                                                       const arma::colvec &arg)
{
    return warpingSet.dissimilarityPointer->GetDistance(
            warpingSet.inputGrid1,
            warpingSet.inputGrid2,
            warpingSet.inputValues1,
            warpingSet.inputValues2
    );
}

#ifndef MEDOIDCENTERCLASS_H
#define MEDOIDCENTERCLASS_H

#include "baseCenterClass.h"

/// Medoid center method finds the real medoid center.
class MedoidCenterMethod : public BaseCenterMethod
{
public:
    CenterType GetCenter(
            const arma::mat& inputGrid,
            const arma::cube& inputValues,
            const std::shared_ptr<BaseDissimilarityFunction>& dissimilarityPointer
    );

    CenterType GetCenter(
            const arma::mat& inputGrid,
            const arma::cube& inputValues,
            const std::shared_ptr<BaseDissimilarityFunction>& dissimilarityPointer,
            unsigned int nbThreads
    );
};

#endif /* MEDOIDCENTERCLASS_H */

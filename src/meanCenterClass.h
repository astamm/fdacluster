#ifndef MEANCENTERCLASS_H
#define MEANCENTERCLASS_H

#include "baseCenterClass.h"

/// Mean center method compute an approximation by lowess.
class MeanCenterMethod : public BaseCenterMethod
{
public:
    CenterType GetCenter(
            const arma::mat& inputGrid,
            const arma::cube& inputValues,
            const std::shared_ptr<Dissimilarity>& dissimilarityPointer
    );
};

#endif /* MEANCENTERCLASS_H */

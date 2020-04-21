#ifndef MEANCENTERCLASS_H
#define MEANCENTERCLASS_H

#include "baseCenterClass.h"

class MeanCenterMethod : public BaseCenterMethod
{
public:
    CenterType GetCenter(
            const arma::mat& inputGrid,
            const arma::cube& inputValues,
            const std::shared_ptr<BaseDissimilarityFunction>& dissimilarityPointer
    );
};

#endif /* MEANCENTERCLASS_H */

#ifndef MEDIANCENTERCLASS_H
#define MEDIANCENTERCLASS_H

#include "baseCenterClass.h"

class MedianCenterMethod : public BaseCenterMethod
{
public:
  CenterType GetCenter(
      const arma::mat& inputGrid,
      const arma::cube& inputValues,
      const std::shared_ptr<BaseDissimilarityFunction>& dissimilarityPointer
  );
};

#endif /* MEDIANCENTERCLASS_H */

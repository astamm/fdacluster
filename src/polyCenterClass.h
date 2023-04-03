#ifndef POLYCENTERCLASS_H
#define POLYCENTERCLASS_H

#include "baseCenterClass.h"

class PolyCenterMethod : public BaseCenterMethod
{
public:
  PolyCenterMethod() {}

  CenterType GetCenter(
      const arma::mat& inputGrid,
      const arma::cube& inputValues,
      const std::shared_ptr<BaseDissimilarityFunction>& dissimilarityPointer
  );
};

#endif /* POLYCENTERCLASS_H */

#ifndef POLYCENTERCLASS_H
#define POLYCENTERCLASS_H

#include "baseCenterClass.h"

class PolyCenterMethod : public BaseCenterMethod
{
public:
  PolyCenterMethod()
  {
    m_PolynomialOrder = 4;
  }

  CenterType GetCenter(
      const arma::mat& inputGrid,
      const arma::cube& inputValues,
      const std::shared_ptr<BaseDissimilarityFunction>& dissimilarityPointer
  );

private:
  unsigned int m_PolynomialOrder;
};

#endif /* POLYCENTERCLASS_H */

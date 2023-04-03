#ifndef LOWESSCENTERCLASS_H
#define LOWESSCENTERCLASS_H

#include "baseCenterClass.h"

class LowessCenterMethod : public BaseCenterMethod
{
public:
  LowessCenterMethod()
  {
    m_StatsPackage = Rcpp::Environment("package:stats");
  }

  CenterType GetCenter(
      const arma::mat& inputGrid,
      const arma::cube& inputValues,
      const std::shared_ptr<BaseDissimilarityFunction>& dissimilarityPointer
  );

private:
  Rcpp::Environment m_StatsPackage;
};

#endif /* LOWESSCENTERCLASS_H */

#ifndef LOWESSCENTERCLASS_H
#define LOWESSCENTERCLASS_H

#include "baseCenterClass.h"

class LowessCenterMethod : public BaseCenterMethod
{
public:
  LowessCenterMethod()
  {
    m_SpanValue = 0.1;
    m_StatsPackage = Rcpp::Environment("package:stats");
  }

  CenterType GetCenter(
      const arma::mat& inputGrid,
      const arma::cube& inputValues,
      const std::shared_ptr<BaseDissimilarityFunction>& dissimilarityPointer
  );

private:
  double m_SpanValue;
  Rcpp::Environment m_StatsPackage;
};

#endif /* LOWESSCENTERCLASS_H */

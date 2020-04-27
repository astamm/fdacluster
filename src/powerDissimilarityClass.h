#ifndef POWERDISSIMILARITYCLASS_H
#define POWERDISSIMILARITYCLASS_H

#include "baseDissimilarityClass.h"

/// Power Distance: properly weighted L2 distance
/// for consistency with power warping function
class PowerDissimilarityFunction : public BaseDissimilarityFunction
{
public:
  PowerDissimilarityFunction()
  {
    m_DomainLowerBound = 0.0;
    m_DomainUpperBound = 1.0;
  }

  void SetDomainLowerBound(const double &val) {m_DomainLowerBound = val;}
  void SetDomainUpperBound(const double &val) {m_DomainUpperBound = val;}

  double GetDistance(
      const arma::rowvec& grid1,
      const arma::rowvec& grid2,
      const arma::mat& values1,
      const arma::mat& values2
  );

private:
  double m_DomainLowerBound;
  double m_DomainUpperBound;
};

#endif /* POWERDISSIMILARITYCLASS_H */

#ifndef POWERWARPINGCLASS_H
#define POWERWARPINGCLASS_H

#include "baseWarpingClass.h"

class PowerWarpingFunction : public BaseWarpingFunction
{
  /**
   * A power transformation of the abscissa x given
   *  a power parameter p is: a + (b-a) * ((s-a) / (b-a))^p.
   */

public:
  PowerWarpingFunction()
  {
    m_DomainLowerBound = 0.0;
    m_DomainUpperBound = 1.0;
  }

  void SetDomainLowerBound(const double &val) {m_DomainLowerBound = val;}
  void SetDomainUpperBound(const double &val) {m_DomainUpperBound = val;}

  unsigned int GetNumberOfParameters();

  arma::rowvec GetInitialPoint();

  arma::mat ApplyWarping(
      const arma::mat &inputGrids,
      const arma::mat &warpingParameters
  );

  void SetParameterBounds(
      const arma::rowvec &warpingOptions,
      const arma::mat &inputGrids
  );

  arma::mat GetFinalWarping(
      const arma::cube &warpingParametersContainer,
      const arma::urowvec &observationMemberships,
      const arma::urowvec &clusterIndices
  );

  void Normalize(
      arma::mat &warpingParameters,
      const arma::urowvec &clusterIndices,
      const arma::urowvec &observationMemberships
  );

private:
  double m_DomainLowerBound;
  double m_DomainUpperBound;
};

#endif /* POWERWARPINGCLASS_H */

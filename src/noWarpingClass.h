#ifndef NOWARPINGCLASS_H
#define NOWARPINGCLASS_H

#include "baseWarpingClass.h"

class NoWarpingFunction : public BaseWarpingFunction
{
  /**
   * The applied transformation is the identity:
   */

public:
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
};

#endif /* NOWARPINGCLASS_H */

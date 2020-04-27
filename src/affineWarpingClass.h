#ifndef AFFINEWARPINGCLASS_H
#define AFFINEWARPINGCLASS_H

#include "baseWarpingClass.h"

class AffineWarpingFunction : public BaseWarpingFunction
{
  /**
   * An affine transformation of the abscissa x given
   *  dilation and shift (d,s) is: d * x + s .
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

#endif /* AFFINEWARPINGCLASS_H */

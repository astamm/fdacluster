#ifndef SHIFTWARPINGCLASS_H
#define SHIFTWARPINGCLASS_H

#include "baseWarpingClass.h"

class ShiftWarpingFunction : public BaseWarpingFunction
{
  /**
   * A shift transformation of the abscissa x given
   * the shift ( s ) is: x + s .
   */

public:
  unsigned int GetNumberOfParameters();

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

  double GetDissimilarityAfterWarping(
      const WarpingSet &warpingSet,
      const arma::rowvec &warpingParameters
  );
};

#endif /* SHIFTWARPINGCLASS_H */

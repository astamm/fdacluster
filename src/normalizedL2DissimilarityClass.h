#ifndef NORMALIZEDL2DISSIMILARITYCLASS_H
#define NORMALIZEDL2DISSIMILARITYCLASS_H

#include "baseDissimilarityClass.h"

/// Normalized L2 Distance
class NormalizedL2DissimilarityFunction : public BaseDissimilarityFunction
{
public:
  double GetDistance(
      const arma::rowvec& grid1,
      const arma::rowvec& grid2,
      const arma::mat& values1,
      const arma::mat& values2
  );
};

#endif /* NORMALIZEDL2DISSIMILARITYCLASS_H */

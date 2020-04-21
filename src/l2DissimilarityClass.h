#ifndef L2DISSIMILARITYCLASS_H
#define L2DISSIMILARITYCLASS_H

#include "baseDissimilarityClass.h"

/// L2 Distance
class L2DissimilarityFunction : public BaseDissimilarityFunction
{
public:
  double GetDistance(
      const arma::rowvec& grid1,
      const arma::rowvec& grid2,
      const arma::mat& values1,
      const arma::mat& values2
  );
};

#endif /* L2DISSIMILARITYCLASS_H */

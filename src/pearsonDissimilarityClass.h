#ifndef PEARSONDISSIMILARITYCLASS_H
#define PEARSONDISSIMILARITYCLASS_H

#include "baseDissimilarityClass.h"

/// Pearson similarity
class PearsonDissimilarityFunction : public BaseDissimilarityFunction
{
  /**
   *The function is returned negative to be a dissimilarity and
   * fit the minimization of the algorithm.
   */

public:
  double GetDistance(
      const arma::rowvec& grid1,
      const arma::rowvec& grid2,
      const arma::mat& values1,
      const arma::mat& values2
  );
};

#endif /* PEARSONDISSIMILARITYCLASS_H */

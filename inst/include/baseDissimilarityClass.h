#ifndef BASEDISSIMILARITYCLASS_H
#define BASEDISSIMILARITYCLASS_H

#include <RcppArmadillo.h>

/// Element with approximated abscissa and approximated functions to compute dissimilarity.
struct FunctionPairType
{
  arma::rowvec Grid;
  arma::mat Values1;
  arma::mat Values2;
};

/// Base class for all the available dissimilarity.
class BaseDissimilarityFunction
{
public:
  BaseDissimilarityFunction() {}
  virtual ~BaseDissimilarityFunction() {}

  /// compute dissimilarity method different for each derived class
  /**
   * @param[grid1] evaluation grid of the first function;
   * @param[grid2] evaluation grid of the second function;
   * @param[values1] values of the first function;
   * @param[values2] values of the second function;
   *
   * @return dissimilarity between the two input functions.
   */
  virtual double GetDistance(
      const arma::rowvec& grid1,
      const arma::rowvec& grid2,
      const arma::mat& values1,
      const arma::mat& values2
  ) = 0;

  /// compute common grid and approximations to compute dissimilarity
  /**
   * @param[grid1] evaluation grid of the first function;
   * @param[grid2] evaluation grid of the second function;
   * @param[values1] values of the first function;
   * @param[values2] values of the second function;
   *
   * @return FunctionPairType object containing the two functions evaluated on their common grid.
   */
  FunctionPairType GetComparableFunctions(
      const arma::rowvec& grid1,
      const arma::rowvec& grid2,
      const arma::mat& values1,
      const arma::mat& values2
  );
};

#endif /* BASEDISSIMILARITYCLASS_H */

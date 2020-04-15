#ifndef BASEWARPINGCLASS_H
#define BASEWARPINGCLASS_H

#include "baseDissimilarityClass.h"

#include <RcppArmadillo.h>

/**
 * Input for warp member function.
 * It contains abscisa and values of the functions to warp and it is returned by
 * the set_function member function.
 */
struct WarpingSet
{
  arma::rowvec inputGrid1;
  arma::rowvec inputGrid2;
  arma::mat inputValues1;
  arma::mat inputValues2;
  std::shared_ptr<BaseDissimilarityFunction> dissimilarityPointer;
};

/// Base class for all warping functions
/**
 * From this class all the warping available are derived.
 */
class BaseWarpingFunction
{
public:
  /// Member to create WarpingSet.
  /**
   * @param[grid1] evaluation grid of the first function;
   * @param[grid2] evaluation grid of the second function;
   * @param[values1] values of the first function on the corresponding evaluation grid;
   * @param[values2] values of the second function on the corresponding evaluation grid;
   * @param[dissimilarity] pointer to dissimilarity object;
   *
   * @return a WarpingSet object.
   */
  WarpingSet SetInputData(
      const arma::rowvec &grid1,
      const arma::rowvec &grid2,
      const arma::mat &values1,
      const arma::mat &values2,
      const std::shared_ptr<BaseDissimilarityFunction> &dissimilarity
  );

  arma::rowvec GetParameterLowerBounds() {return m_ParameterLowerBounds;}
  arma::rowvec GetParameterUpperBounds() {return m_ParameterUpperBounds;}

  /// Apply warping to a matrix.
  /**
   * @param[x] abscissa to warp;
   * @param[par] warping parameters to apply;
   *
   * return warped abscissas.
   */
  virtual arma::mat ApplyWarping(const arma::mat &x, const arma::mat &par) = 0;

  /// Return number of parameters.
  virtual unsigned int GetNumberOfParameters() = 0;

  /// Set bounds given the input option different for each warping function.
  /**
   * @param[warpingOptions] input warping option;
   * @param[x] abscissa to warp.
   */
  virtual void SetParameterBounds(const arma::rowvec &warpingOptions, const arma::mat &x) = 0;

  /// Set bounds given in a matrix.
  /**
   * @param[bounds] bounds already computed;
   */
  void SetParameterBounds(const arma::mat &bounds);

  /// Compute final warping.
  /**
   * @param[parameters_vec] warping's parameters of each iteration;
   * @param[labels] final labels;
   * @param[ict] index current clusters;
   *
   * @return a matrix with the total warping parameters applied.
   */
  virtual arma::mat GetFinalWarping(
      const arma::cube &parameters_vec,
      const arma::urowvec &labels,
      const arma::urowvec &ict
  ) = 0;

  /// Normalize the warping parameters computed by clusters.
  /**
   * @param[par] warping parameters computed;
   * @param[ict] index current clusters;
   * @param[labels] current labels;
   */
  virtual void Normalize(
      arma::mat &par,
      const arma::urowvec &ict,
      const arma::urowvec &labels
  ) = 0;

  /// Compute dissimilarity after warp for optimization.
  /**
   * @param[w_set] warping_set element with the functions to warp.
   * @param[arg] best parameters that will be computed.
   */
  virtual double GetDissimilarityAfterWarping(
      const WarpingSet &warpingSet,
      const arma::colvec &arg) = 0;

protected:
  arma::rowvec m_ParameterLowerBounds;
  arma::rowvec m_ParameterUpperBounds;
};

#endif /* BASEWARPINGCLASS_H */

#ifndef BASEWARPINGCLASS_H
#define BASEWARPINGCLASS_H

#include "baseDissimilarityClass.h"

#include <RcppArmadillo.h>
#include <memory>

/**
 * Input for warp member function.
 * It contains abscissa and values of the functions to warp and it is returned by
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
  BaseWarpingFunction()
  {
    m_ParameterLowerBounds.reset();
    m_ParameterUpperBounds.reset();
  }

  virtual ~BaseWarpingFunction() {}

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

  /// Return number of parameters.
  virtual unsigned int GetNumberOfParameters() = 0;

  /// Return a initial guess at the parameters to start optimization with
  virtual arma::rowvec GetInitialPoint() = 0;

  /// Apply warping to a matrix.
  /**
   * @param[inputGrids] abscissa to warp;
   * @param[warpingParameters] warping parameters to apply;
   *
   * return warped abscissas.
   */
  virtual arma::mat ApplyWarping(
      const arma::mat &inputGrids,
      const arma::mat &warpingParameters
  ) = 0;

  /// Set bounds given the input option different for each warping function.
  /**
   * @param[warpingOptions] input warping option;
   * @param[inputGrids] abscissa to warp.
   */
  virtual void SetParameterBounds(
      const arma::rowvec &warpingOptions,
      const arma::mat &inputGrids
  ) = 0;

  /// Set bounds given in a matrix.
  /**
   * @param[bounds] bounds already computed;
   */
  void SetParameterBounds(const arma::mat &bounds);

  /// Compute final warping.
  /**
   * @param[warpingParametersContainer] warping parameters of each iteration;
   * @param[observationMemberships] final labels;
   * @param[clusterIndices] index current clusters;
   *
   * @return a matrix with the total warping parameters applied.
   */
  virtual arma::mat GetFinalWarping(
      const arma::cube &warpingParametersContainer,
      const arma::urowvec &observationMemberships,
      const arma::urowvec &clusterIndices
  ) = 0;

  /// Normalize the warping parameters computed by clusters.
  /**
   * @param[warpingParameters] warping parameters computed;
   * @param[clusterIndices] index current clusters;
   * @param[observationMemberships] current labels;
   */
  virtual void Normalize(
      arma::mat &warpingParameters,
      const arma::urowvec &clusterIndices,
      const arma::urowvec &observationMemberships
  ) = 0;

  /// Compute dissimilarity after warp for optimization.
  /**
   * @param[warpingSet] warping_set element with the functions to warp.
   * @param[warpingParameters] best parameters that will be computed.
   */
  double GetDissimilarityAfterWarping(
      const WarpingSet &warpingSet,
      const arma::rowvec &warpingParameters
  );

protected:
  arma::rowvec m_ParameterLowerBounds;
  arma::rowvec m_ParameterUpperBounds;
};

#endif /* BASEWARPINGCLASS_H */

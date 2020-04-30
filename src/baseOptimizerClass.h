#ifndef BASEOPTIMIZERCLASS_H
#define BASEOPTIMIZERCLASS_H

#include "baseWarpingClass.h"

#include <RcppArmadillo.h>
#include <nloptrAPI.h>

class BaseOptimizerFunction
{
public:
  struct AlignToTemplateData
  {
    arma::rowvec inputGrid;
    arma::mat inputValues;
    arma::rowvec templateGrid;
    arma::mat templateValues;
    std::shared_ptr<BaseDissimilarityFunction> dissimilarityPointer;
    std::shared_ptr<BaseWarpingFunction> warpingPointer;
    double penalizationWeight;
  };

  struct CenterTemplateData
  {
    arma::rowvec templateGrid;
    arma::mat templateValues;
    arma::mat inputGrids;
    arma::cube inputValues;
    std::shared_ptr<BaseDissimilarityFunction> dissimilarityPointer;
    std::shared_ptr<BaseWarpingFunction> warpingPointer;
    double penalizationWeight;
  };

  BaseOptimizerFunction()
  {
    m_ParameterRelativeTolerance = 1.0e-4;
    m_PenalizationWeight = 0.0;
  }

  virtual ~BaseOptimizerFunction() {}

  void SetPenalizationWeight(const double &val) {m_PenalizationWeight = val;}

  virtual nlopt_opt GetOptimizer(const unsigned int numberOfParameters) = 0;

  double AlignToTemplate(
      arma::rowvec &initialParameters,
      const arma::rowvec &inputGrid,
      const arma::mat &inputValues,
      const arma::rowvec &templateGrid,
      const arma::mat &templateValues,
      const std::shared_ptr<BaseDissimilarityFunction> &dissimilarityPointer,
      const std::shared_ptr<BaseWarpingFunction> &warpingPointer
  );

  double CenterTemplate(
      arma::rowvec &initialParameters,
      const arma::rowvec &templateGrid,
      const arma::mat &templateValues,
      const arma::mat &inputGrids,
      const arma::cube &inputValues,
      const std::shared_ptr<BaseDissimilarityFunction> &dissimilarityPointer,
      const std::shared_ptr<BaseWarpingFunction> &warpingPointer
  );

protected:
  static double AlignToTemplateCostFunction(
      unsigned n,
      const double *x,
      double *grad,
      void *data
  );
  static double CenterTemplateCostFunction(
      unsigned n,
      const double *x,
      double *grad,
      void *data
  );

private:
  double m_ParameterRelativeTolerance;
  double m_PenalizationWeight;
};

#endif /* BASEOPTIMIZERCLASS_H */

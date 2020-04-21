#ifndef BASEOPTIMIZERCLASS_H
#define BASEOPTIMIZERCLASS_H

#include "baseWarpingClass.h"

#include <RcppArmadillo.h>
#include <nloptrAPI.h>

class BaseOptimizerFunction
{
public:
  struct CostFunctionData
  {
    std::shared_ptr<BaseWarpingFunction> warpingPointer;
    WarpingSet warpingSet;
  };

  BaseOptimizerFunction()
  {
    m_ParameterRelativeTolerance = 1.0e-4;
  }

  virtual ~BaseOptimizerFunction() {}

  virtual nlopt_opt GetOptimizer(const unsigned int numberOfParameters) = 0;

  double Optimize(
      arma::rowvec &initialParameters,
      const std::shared_ptr<BaseWarpingFunction> &warpingPointer,
      const WarpingSet &warpingSet
  );

protected:
  static double GetValue(unsigned n, const double *x, double *grad, void *data);

private:
  double m_ParameterRelativeTolerance;
};

#endif /* BASEOPTIMIZERCLASS_H */

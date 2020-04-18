#ifndef BASEOPTIMIZERCLASS_H
#define BASEOPTIMIZERCLASS_H

#include "baseWarpingClass.h"

#include <RcppArmadillo.h>

struct CostFunctionData
{
  std::shared_ptr<BaseWarpingFunction> warpingPointer;
  WarpingSet warpingSet;
};

double GetCostFunctionValue(unsigned n, const double *x, double *grad, void *data);

class BaseOptimizerFunction
{
public:
  BaseOptimizerFunction() {}
  virtual ~BaseOptimizerFunction() {}

  virtual double Optimize(
      arma::rowvec &initialParameters,
      const std::shared_ptr<BaseWarpingFunction> &warpingPointer,
      const WarpingSet &warpingSet
  ) = 0;
};

#endif /* BASEOPTIMIZERCLASS_H */

#ifndef BASEOPTIMIZERCLASS_H
#define BASEOPTIMIZERCLASS_H

#include "baseWarpingClass.h"

#include <RcppArmadillo.h>
#include <dlib/optimization/find_optimal_parameters.h>

class BaseOptimizerFunction
{
public:
  using ParametersType = dlib::matrix<double,0,1>;

  virtual double Optimize(
      arma::rowvec &initialParameters,
      const std::shared_ptr<BaseWarpingFunction> &warpingFunction,
      const std::function<double(arma::rowvec)> &costFunction
  ) = 0;
};

#endif /* BASEOPTIMIZERCLASS_H */

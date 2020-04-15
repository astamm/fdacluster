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
            arma::rowvec &arg,
            std::shared_ptr<BaseWarpingFunction> &warpingFunction,
            std::function<double(arma::colvec)> fun) = 0;
};

#endif /* BASEOPTIMIZERCLASS_H */

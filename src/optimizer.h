#ifndef MY_OPTIMIZER_HPP
#define MY_OPTIMIZER_HPP

#include <RcppArmadillo.h>
#include "warping.h"

//
//  optimizer
//

class OptimizerMethod{
  public:
    virtual double optimize(colvec& arg, std::shared_ptr<WarpingFunction>& pfunc, std::function<double(colvec)> fun)=0;
};


class Bobyqa: public OptimizerMethod{
public:
  virtual double optimize(colvec& arg, std::shared_ptr<WarpingFunction>& pfunc, std::function<double(colvec)> fun);
};

#endif

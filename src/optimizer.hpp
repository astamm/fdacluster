#ifndef MY_OPTIMIZER_HPP
#define MY_OPTIMIZER_HPP

#include <RcppArmadillo.h>
#include "warping.hpp"

//
//  optimizer
//

class OptimizerMethod{
  public:
    virtual double optimize(colvec& arg, std::shared_ptr<WarpingFunction>& pfunc, std::function<double(colvec)> fun)=0;
};



class Bobyca: public OptimizerMethod{
public:
  virtual double optimize(colvec& arg, std::shared_ptr<WarpingFunction>& pfunc, std::function<double(colvec)> fun);
};

#endif

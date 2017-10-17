#ifndef MY_OPTIMIZER_HPP
#define MY_OPTIMIZER_HPP

#include <RcppArmadillo.h>
#include "warping.h"

//
//  optimizer
//

/// Optimizer Base class
class OptimizerMethod{
  public:
    ///optimer member function
    /**
     * @param[arg] warping's warping to be optimized returned by reference;
     * @param[pfunc] pointer to warping object;
     * @param[fun] functor(colvec) to be optimized, input customizable.
     *
     * return optimal dissimilarity.
     */
    virtual double optimize(colvec& arg, std::shared_ptr<WarpingFunction>& pfunc, std::function<double(colvec)> fun)=0;
};


class Bobyqa: public OptimizerMethod{
public:
  virtual double optimize(colvec& arg, std::shared_ptr<WarpingFunction>& pfunc, std::function<double(colvec)> fun);
};

#endif

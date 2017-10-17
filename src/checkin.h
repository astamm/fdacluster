#ifndef CHECK_IN_HPP
#define CHECK_IN_HPP

#include<RcppArmadillo.h>
#include "dissimilarity.h"
#include "warping.h"
#include "center_methods.h"
#include "optimizer.h"

bool checkIn(const arma::mat& x,
             const arma::cube& y,
             const std::string warping_method,
             const std::string center_method,
             const std::string similarity_method,
             const std::string optim_method,
             const util::SharedFactory<WarpingFunction> a,
             const util::SharedFactory<Dissimilarity> b,
             const util::SharedFactory<CenterMethod> c,
             const util::SharedFactory<OptimizerMethod> d,
             const arma::rowvec& warping_opt,
             const arma::rowvec& center_opt,
             const arma::urowvec par_opt);
#endif

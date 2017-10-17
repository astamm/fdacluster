#ifndef CHECK_IN_HPP
#define CHECK_IN_HPP

#include<RcppArmadillo.h>
#include "dissimilarity.h"
#include "warping.h"
#include "center_methods.h"
#include "optimizer.h"

void checkIn(const arma::mat& x,
             const arma::cube& y,
             const std::string warping_method,
             const std::string center_method,
             const std::string similarity_method,
             const std::string optim_method,
             const util::SharedFactory<WarpingFunction> warfac,
             const util::SharedFactory<Dissimilarity> disfac,
             const util::SharedFactory<CenterMethod> cenfac,
             const util::SharedFactory<OptimizerMethod> optfac,
             const arma::rowvec& warping_opt,
             const arma::rowvec& center_opt,
             const arma::urowvec par_opt,
             bool show_iter);
#endif

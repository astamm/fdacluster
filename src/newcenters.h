#ifndef NEW_CENTERS_HPP
#define NEW_CENTERS_HPP

#include <RcppArmadillo.h>

#include "dissimilarity.h"
#include "center_methods.h"

/// Compute new centers
/**
 * It's a function that handle the computation of new centers.
 */

void newCenters(const arma::mat& x_reg,
                      const arma::cube& y,
                      const arma::rowvec& x_out,
                      std::shared_ptr<Dissimilarity>& dissim,
                      std::shared_ptr<CenterMethod>& cen,
                      const arma::urowvec& par_opt,
                      cube& templates,
                      const urowvec& ict,
                      const urowvec& labels,
                      const bool show_iter
                     );

#endif

#ifndef NEWCENTERSMETHOD_H
#define NEWCENTERSMETHOD_H

#include "baseDissimilarityClass.h"
#include "baseCenterClass.h"

/// Compute new centers
/**
 * It's a function that handle the computation of new centers.
 */

void newCenters(const arma::mat& x_reg,
                const arma::cube& y,
                const arma::rowvec& x_out,
                std::shared_ptr<BaseDissimilarityFunction>& dissim,
                std::shared_ptr<BaseCenterMethod>& cen,
                const arma::urowvec& par_opt,
                arma::cube& templates,
                const arma::urowvec& ict,
                const arma::urowvec& labels,
                const bool show_iter
);

#endif /* NEWCENTERSMETHOD_H */

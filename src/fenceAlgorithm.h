#ifndef FENCEALGORITHM_H
#define FENCEALGORITHM_H

#include "baseWarpingClass.h"
#include "baseOptimizerClass.h"

/// Remove warping outliers
/** It's an optional check that can be activated by input. After each computation
 *  of best warping if the computed parameters are outliers they are recomputed with
 *  stricter bounds. It\' s computational less expensive way to have lighter warping
 *  insted of decrease input bound (for example: max_shift and max_dilation).
 */
void iterativeFence(
        arma::mat parameters,
        const arma::uword iter,
        arma::urowvec& labels,
        arma::rowvec& index,
        std::shared_ptr<BaseWarpingFunction>& warping,
        std::shared_ptr<BaseOptimizerFunction>& optimizer,
        const arma::cube& templates,
        const arma::mat& x_reg,
        const arma::cube& y,
        const arma::rowvec& x_out,
        std::shared_ptr<BaseDissimilarityFunction>& dissim,
        const arma::urowvec& ict,
        const bool show_iter);

#endif /* FENCEALGORITHM_H */

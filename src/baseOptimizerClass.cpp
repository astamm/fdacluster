#include "baseOptimizerClass.h"

double BaseOptimizerFunction::GetValue(unsigned n, const double *x, double *grad, void *data)
{
  CostFunctionData *d = (CostFunctionData *) data;

  arma::rowvec params(n);
  for (unsigned int i = 0;i < n;++i)
    params(i) = x[i];

  return d->warpingPointer->GetDissimilarityAfterWarping(d->warpingSet, params);
}

double BaseOptimizerFunction::Optimize(arma::rowvec &initialParameters,
                                       const std::shared_ptr<BaseWarpingFunction> &warpingPointer,
                                       const WarpingSet &warpingSet)
{
    unsigned int numberOfParameters = warpingPointer->GetNumberOfParameters();
    nlopt_opt optimizer = this->GetOptimizer(numberOfParameters);

    arma::rowvec lowerBounds = warpingPointer->GetParameterLowerBounds();
    arma::rowvec upperBounds = warpingPointer->GetParameterUpperBounds();
    initialParameters = warpingPointer->GetInitialPoint();

    CostFunctionData data;
    data.warpingPointer = warpingPointer;
    data.warpingSet = warpingSet;

    if (initialParameters.size() == 0)
        return this->GetValue(numberOfParameters, &(initialParameters[0]), NULL, &data);

    nlopt_set_lower_bounds(optimizer, &(lowerBounds[0]));
    nlopt_set_upper_bounds(optimizer, &(upperBounds[0]));

    nlopt_set_min_objective(optimizer, this->GetValue, &data);
    nlopt_set_xtol_rel(optimizer, m_ParameterRelativeTolerance);

    double fVal;
    int exitCode = nlopt_optimize(optimizer, &(initialParameters[0]), &fVal);

    nlopt_destroy(optimizer);

    if (exitCode < 0)
        Rcpp::stop("NLOPT optimization failed.");

    return fVal;
}

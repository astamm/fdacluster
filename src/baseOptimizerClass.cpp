#include "baseOptimizerClass.h"

double GetCostFunctionValue(unsigned n, const double *x, double *grad, void *data)
{
  CostFunctionData *d = (CostFunctionData *) data;

  arma::rowvec params(n);
  for (unsigned int i = 0;i < n;++i)
    params(i) = x[i];

  return d->warpingPointer->GetDissimilarityAfterWarping(d->warpingSet, params);
}

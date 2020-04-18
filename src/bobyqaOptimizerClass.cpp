#include "bobyqaOptimizerClass.h"

#include <nloptrAPI.h>

double BobyqaOptimizerFunction::Optimize(arma::rowvec &initialParameters,
                                         const std::shared_ptr<BaseWarpingFunction> &warpingPointer,
                                         const WarpingSet &warpingSet)
{
    unsigned int numberOfParameters = warpingPointer->GetNumberOfParameters();
    nlopt_opt optimizer = nlopt_create(NLOPT_LN_BOBYQA, numberOfParameters);

    arma::rowvec lowerBounds = warpingPointer->GetParameterLowerBounds();
    arma::rowvec upperBounds = warpingPointer->GetParameterUpperBounds();
    initialParameters = (lowerBounds + upperBounds) / 2.0;

    CostFunctionData data;
    data.warpingPointer = warpingPointer;
    data.warpingSet = warpingSet;

    if (initialParameters.size() == 0)
        return GetCostFunctionValue(numberOfParameters, &(initialParameters[0]), NULL, &data);

    nlopt_set_lower_bounds(optimizer, &(lowerBounds[0]));
    nlopt_set_upper_bounds(optimizer, &(upperBounds[0]));

    nlopt_set_min_objective(optimizer, GetCostFunctionValue, &data);
    nlopt_set_xtol_rel(optimizer, 1e-4);

    double fVal;
    int exitCode = nlopt_optimize(optimizer, &(initialParameters[0]), &fVal);

    nlopt_destroy(optimizer);

    if (exitCode < 0)
        Rcpp::stop("NLOPT optimization failed.");

    return fVal;
}

#include "bobyqaOptimizerClass.h"

double BobyqaOptimizerFunction::Optimize(arma::rowvec &initialParameters,
                                         const std::shared_ptr<BaseWarpingFunction> &warpingPointer,
                                         const WarpingSet &warpingSet)
{
    unsigned int numberOfParameters = warpingPointer->GetNumberOfParameters();
    arma::rowvec lowerBounds = warpingPointer->GetParameterLowerBounds();
    arma::rowvec upperBounds = warpingPointer->GetParameterUpperBounds();
    initialParameters = (lowerBounds + upperBounds) / 2.0;
    arma::rowvec d = upperBounds - lowerBounds;

    ParametersType dlibInitialParameters(numberOfParameters);
    ParametersType dlibLowerBounds(numberOfParameters);
    ParametersType dlibUpperBounds(numberOfParameters);

    for (unsigned int i = 0;i < numberOfParameters;++i)
    {
        dlibInitialParameters(i) = initialParameters(i);
        dlibLowerBounds(i) = lowerBounds(i);
        dlibUpperBounds(i) = upperBounds(i);
    }

    double radius = d.min() / 2.0 - m_EpsilonValue;

    auto dlibCostFunction = [&warpingPointer, &warpingSet] (const ParametersType& dlibParams)
    {
        unsigned int numberOfParameters = dlibParams.nr();
        arma::rowvec params(numberOfParameters);
        for (unsigned int i = 0;i < numberOfParameters;++i)
            params(i) = dlibParams(i);

        return warpingPointer->GetDissimilarityAfterWarping(warpingSet, params);
    };

    if (initialParameters.size() == 0)
        return dlibCostFunction(dlibInitialParameters);

    find_optimal_parameters(
        radius,
        m_EpsilonValue,
        100,
        dlibInitialParameters,
        dlibLowerBounds,
        dlibUpperBounds,
        dlibCostFunction
    );

    for (unsigned int i = 0;i < numberOfParameters;++i)
        initialParameters(i) = dlibInitialParameters(i);

    return dlibCostFunction(dlibInitialParameters);
}



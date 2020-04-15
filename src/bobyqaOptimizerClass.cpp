#include "bobyqaOptimizerClass.h"

double BobyqaOptimizerFunction::Optimize(arma::rowvec &initialParameters,
                                         const std::shared_ptr<BaseWarpingFunction> &warpingFunction,
                                         const std::function<double(arma::rowvec)> &costFunction)
{
    unsigned int numberOfParameters = warpingFunction->GetNumberOfParameters();
    arma::rowvec lowerBounds = warpingFunction->GetParameterLowerBounds();
    arma::rowvec upperBounds = warpingFunction->GetParameterUpperBounds();
    initialParameters = (lowerBounds + upperBounds) / 2.0;

    if (initialParameters.size() == 0)
        return costFunction(initialParameters);

    auto dlibCostFunction = [&costFunction] (const ParametersType& dlibParams)
    {
        unsigned int numberOfParameters = dlibParams.nr();

        arma::rowvec params(numberOfParameters);
        for (unsigned int i = 0;i < numberOfParameters;++i)
            params(i) = dlibParams(i);

        return costFunction(params);
    };

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

    return costFunction(initialParameters);
}



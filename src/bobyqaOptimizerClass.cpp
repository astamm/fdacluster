#include "bobyqaOptimizerClass.h"

double BobyqaOptimizerFunction::Optimize(arma::colvec &arg,
                                         std::shared_ptr<BaseWarpingFunction> &warpingFunction,
                                         std::function<double(arma::colvec)> fun)
{
    auto dlibCostFunction = [&fun] (const ParametersType& argt)
    {
        unsigned int numberOfParameters = argt.nr();

        arma::colvec args(numberOfParameters);
        for (unsigned int i = 0;i < numberOfParameters;++i)
            args(i) = argt(i);

        return fun(args);
    };

    unsigned int numberOfParameters = warpingFunction->GetNumberOfParameters();
    arma::rowvec lowerBounds = warpingFunction->GetParameterLowerBounds();
    arma::rowvec upperBounds = warpingFunction->GetParameterUpperBounds();
    arma::rowvec s = (lowerBounds + upperBounds) / 2.0;
    arma::rowvec d = upperBounds - lowerBounds;

    ParametersType dlibStartingPoint(numberOfParameters);
    ParametersType dlibLowerBounds(numberOfParameters);
    ParametersType dlibUpperBounds(numberOfParameters);

    for (unsigned int i = 0;i < numberOfParameters;++i)
    {
        dlibStartingPoint(i) = s(i);
        dlibLowerBounds(i) = lowerBounds(i);
        dlibUpperBounds(i) = upperBounds(i);
    }

    // gestione caso no aligment numero di parametri da stimare 0
    if (arg.size() == 0)
        return dlibCostFunction(dlibStartingPoint);

    double radius = d.min() / 2.0 - m_EpsilonValue;

    find_optimal_parameters(
        radius,
        m_EpsilonValue,
        100,
        dlibStartingPoint,
        dlibLowerBounds,
        dlibUpperBounds,
        dlibCostFunction
    );

    for (unsigned int i = 0;i < numberOfParameters;++i)
        arg(i) = dlibStartingPoint(i);

    return dlibCostFunction(dlibStartingPoint);
}



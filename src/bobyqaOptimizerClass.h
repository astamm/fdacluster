#ifndef BOBYQAOPTIMIZERCLASS_H
#define BOBYQAOPTIMIZERCLASS_H

#include "baseOptimizerClass.h"

class BobyqaOptimizerFunction : public BaseOptimizerFunction
{
public:
    BobyqaOptimizerFunction()
    {
        m_EpsilonValue = 1.0e-8;
    }

    void SetEpsilonValue(const double &val) {m_EpsilonValue = val;}

    double Optimize(
            arma::rowvec &initialParameters,
            const std::shared_ptr<BaseWarpingFunction> &warpingPointer,
            const WarpingSet &warpingSet
    );

private:
    double m_EpsilonValue;
};

#endif /* BOBYQAOPTIMIZERCLASS_H */

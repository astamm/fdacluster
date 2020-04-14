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
            arma::colvec &arg,
            std::shared_ptr<BaseWarpingFunction> &warpingFunction,
            std::function<double(arma::colvec)> fun);

private:
    double m_EpsilonValue;
};

#endif /* BOBYQAOPTIMIZERCLASS_H */

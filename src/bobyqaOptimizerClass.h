#ifndef BOBYQAOPTIMIZERCLASS_H
#define BOBYQAOPTIMIZERCLASS_H

#include "baseOptimizerClass.h"

class BobyqaOptimizerFunction : public BaseOptimizerFunction
{
public:
    BobyqaOptimizerFunction() {}
    virtual ~BobyqaOptimizerFunction() {}

    nlopt_opt GetOptimizer(const unsigned int numberOfParameters);
};

#endif /* BOBYQAOPTIMIZERCLASS_H */

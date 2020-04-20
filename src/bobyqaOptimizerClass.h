#ifndef BOBYQAOPTIMIZERCLASS_H
#define BOBYQAOPTIMIZERCLASS_H

#include "baseOptimizerClass.h"

class BobyqaOptimizerFunction : public BaseOptimizerFunction
{
public:
    nlopt_opt GetOptimizer(const unsigned int numberOfParameters);
};

#endif /* BOBYQAOPTIMIZERCLASS_H */

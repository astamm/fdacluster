#include "bobyqaOptimizerClass.h"

nlopt_opt BobyqaOptimizerFunction::GetOptimizer(const unsigned int numberOfParameters)
{
    return nlopt_create(NLOPT_LN_BOBYQA, numberOfParameters);
}

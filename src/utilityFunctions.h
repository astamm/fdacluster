#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include <RcppArmadillo.h>
#include <functional> // for std::function
#include <memory> // for std::shared_ptr...
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

/// R table function in C++
/**
 * @param[inputLabels] Vector of input labels to count.
 *
 * @return Table containing the label counts.
 */
std::map<unsigned int, unsigned int> tableCpp(const arma::urowvec &inputLabels);

/// Extract several observations from a cube.
/**
 *  @param[inputData] Data array in arma::cube format.
 *  @param[observationIndices] Indices of the observations to extract.
 *  @return Extracted observations.
 */
arma::cube GetObservations(const arma::cube& inputData, arma::uvec& observationIndices);

/// Factory class
template <class BaseObjectType>
class SharedFactory
{
public:
    using SharedPointerType = std::shared_ptr<BaseObjectType>;
    using RegistryMap = std::unordered_map<std::string, std::function<SharedPointerType()> >;

    // Use this to instantiate the proper Derived class
    SharedPointerType Instantiate(const std::string &name);

    template <class DerivedObjectType> void Register(const std::string &name);

private:
    RegistryMap m_Map;
};

#endif /* UTILITYFUNCTIONS_H */

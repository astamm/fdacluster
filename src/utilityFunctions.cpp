#include "utilityFunctions.h"

//
// tableCpp
//
std::map<unsigned int, unsigned int> tableCpp(const arma::urowvec &inputLabels)
{
    std::map<unsigned int, unsigned int> outputCounts;
    unsigned int numberOfObservations = inputLabels.size();

    for (unsigned int i = 0;i < numberOfObservations;++i)
        ++outputCounts[inputLabels[i]];

    return outputCounts;
}

//
// GetObservations
//
arma::cube GetObservations(const arma::cube& inputData, arma::uvec& observationIndices)
{
    arma::cube outputCube(observationIndices.size(), inputData.n_cols, inputData.n_slices);

    for (unsigned int i = 0;i < observationIndices.size();++i)
        outputCube.row(i) = inputData.row(observationIndices(i));

    return outputCube;
}

template <class BaseObjectType>
typename SharedFactory<BaseObjectType>::SharedPointerType
SharedFactory<BaseObjectType>::Instantiate(const std::string &name)
{
    auto it = m_Map.find(name);
    return it == m_Map.end() ? nullptr : (it->second)();
}

template <class BaseObjectType>
template <class DerivedObjectType>
void SharedFactory<BaseObjectType>::Register(const std::string &name)
{
    m_Map[name] = []()
    {
        return std::make_shared<DerivedObjectType>();
    };
}

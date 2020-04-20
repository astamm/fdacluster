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

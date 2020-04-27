#include "utilityFunctions.h"

std::map<unsigned int, unsigned int> tableCpp(const arma::urowvec &inputLabels)
{
    std::map<unsigned int, unsigned int> outputCounts;
    unsigned int numberOfObservations = inputLabels.size();

    for (unsigned int i = 0;i < numberOfObservations;++i)
        ++outputCounts[inputLabels[i]];

    return outputCounts;
}

arma::cube GetObservations(const arma::cube& inputData, arma::uvec& observationIndices)
{
    arma::cube outputCube(observationIndices.size(), inputData.n_cols, inputData.n_slices);

    for (unsigned int i = 0;i < observationIndices.size();++i)
        outputCube.row(i) = inputData.row(observationIndices(i));

    return outputCube;
}

Rcpp::NumericVector FormatVectorForOutput(const arma::rowvec &inputVector)
{
    Rcpp::NumericVector outputVector = Rcpp::wrap(inputVector);
    outputVector.attr("dim") = R_NilValue;
    return outputVector;
}

Rcpp::NumericVector FormatVectorForOutput(const arma::urowvec &inputVector)
{
    Rcpp::NumericVector outputVector = Rcpp::wrap(inputVector);
    outputVector.attr("dim") = R_NilValue;
    return outputVector;
}

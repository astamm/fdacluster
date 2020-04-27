#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include <RcppArmadillo.h>
#include <map>

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

/// Convert arma::rowvec to a type that converts to vector in R
Rcpp::NumericVector FormatVectorForOutput(const arma::rowvec &inputVector);
Rcpp::NumericVector FormatVectorForOutput(const arma::urowvec &inputVector);

#endif /* UTILITYFUNCTIONS_H */

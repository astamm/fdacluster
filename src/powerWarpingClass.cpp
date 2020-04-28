#include "powerWarpingClass.h"

unsigned int PowerWarpingFunction::GetNumberOfParameters()
{
    return 1;
}

arma::rowvec PowerWarpingFunction::GetInitialPoint()
{
    return { 1 };
}

arma::mat PowerWarpingFunction::ApplyWarping(const arma::mat &inputGrids,
                                             const arma::mat &warpingParameters)
{
    arma::mat outputGrids(inputGrids.n_rows, inputGrids.n_cols);

    for (unsigned int i = 0;i < inputGrids.n_rows;++i)
        outputGrids.row(i) = m_DomainLowerBound + (m_DomainUpperBound - m_DomainLowerBound) * arma::pow((inputGrids.row(i) - m_DomainLowerBound) / (m_DomainUpperBound - m_DomainLowerBound), warpingParameters(i, 0));

    return outputGrids;
}

void PowerWarpingFunction::SetParameterBounds(const arma::rowvec &warpingOptions,
                                              const arma::mat &inputGrids)
{
    double lowerBound = warpingOptions(0);
    if (lowerBound <= 0 || lowerBound >= 1)
        Rcpp::stop("The first entry of the warping option for power warping should be in (0,1) as it corresponds to the power lower bound.");

    double upperBound = warpingOptions(1);
    if (upperBound <= 1)
        Rcpp::stop("The second entry of the warping option for power warping should be in (1,inf) as it corresponds to the power upper bound.");

    m_ParameterLowerBounds = { lowerBound };
    m_ParameterUpperBounds = { upperBound };
}

arma::mat PowerWarpingFunction::GetFinalWarping(const arma::cube &warpingParametersContainer,
                                                const arma::urowvec &observationMemberships,
                                                const arma::urowvec &clusterIndices)
{
    unsigned int numberOfObservations = warpingParametersContainer.n_rows;
    unsigned int numberOfParameters = warpingParametersContainer.n_cols;
    unsigned int numberOfIterations = warpingParametersContainer.n_slices;
    arma::mat warpingParameters(numberOfObservations, numberOfParameters, arma::fill::ones);

    for (unsigned int i = 0;i < numberOfIterations;++i)
        warpingParameters.col(0) %= warpingParametersContainer.slice(i).col(0);

    this->Normalize(warpingParameters, clusterIndices, observationMemberships);

    return warpingParameters;
}

void PowerWarpingFunction::Normalize(arma::mat &warpingParameters,
                                     const arma::urowvec &clusterIndices,
                                     const arma::urowvec &observationMemberships)
{
    arma::uvec observationIndices;
    arma::rowvec meanParameters;
    arma::mat clusterValues;

    for (unsigned int i = 0;i < clusterIndices.size();++i)
    {
        observationIndices = arma::find(observationMemberships == clusterIndices(i));
        clusterValues = warpingParameters.rows(observationIndices);

        double meanPower = std::exp(arma::mean(arma::log(clusterValues.col(0))));
        clusterValues.col(0) /= meanPower;
        warpingParameters.rows(observationIndices) = clusterValues;
    }
}

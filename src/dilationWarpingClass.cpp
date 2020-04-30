#include "dilationWarpingClass.h"

unsigned int DilationWarpingFunction::GetNumberOfParameters()
{
    return 1;
}

arma::rowvec DilationWarpingFunction::GetInitialPoint()
{
    return { 1 };
}

arma::mat DilationWarpingFunction::ApplyWarping(const arma::mat &inputGrids,
                                                const arma::mat &warpingParameters)
{
    arma::mat outputGrids(inputGrids.n_rows, inputGrids.n_cols);

    for (unsigned int i = 0;i < inputGrids.n_rows;++i)
        outputGrids.row(i) = inputGrids.row(i) * warpingParameters(i, 0);

    return outputGrids;
}

void DilationWarpingFunction::SetParameterBounds(const arma::rowvec &warpingOptions,
                                                 const arma::mat &inputGrids)
{
    double dl = warpingOptions(0);
    if (dl < 0 || dl > 1)
        Rcpp::stop("The warping option dl for the dilation parameter should be in [0, 1], as the bounds are computed as [1-dl, 1+dl] centered around the unit dilation.");

    m_ParameterLowerBounds = { 1 - dl };
    m_ParameterUpperBounds = { 1 + dl };
}

arma::mat DilationWarpingFunction::GetFinalWarping(const arma::cube &warpingParametersContainer,
                                                   const arma::urowvec &observationMemberships,
                                                   const arma::urowvec &clusterIndices)
{
    unsigned int numberOfObservations = warpingParametersContainer.n_rows;
    unsigned int numberOfParameters = warpingParametersContainer.n_cols;
    unsigned int numberOfIterations = warpingParametersContainer.n_slices;
    arma::mat warpingParameters(numberOfObservations, numberOfParameters, arma::fill::ones);
    arma::colvec dilationParameters;

    for (unsigned int i = 0;i < numberOfIterations;++i)
    {
        dilationParameters = warpingParametersContainer.slice(i).col(0);
        warpingParameters.col(0) %= dilationParameters;
    }

    this->Normalize(warpingParameters, clusterIndices, observationMemberships);

    return warpingParameters;
}

void DilationWarpingFunction::Normalize(arma::mat &warpingParameters,
                                        const arma::urowvec &clusterIndices,
                                        const arma::urowvec &observationMemberships)
{
    return;

    arma::uvec observationIndices;
    arma::rowvec meanParameters;

    for (unsigned int i = 0;i < clusterIndices.size();++i)
    {
        observationIndices = arma::find(observationMemberships == clusterIndices(i));
        meanParameters = arma::mean(warpingParameters.rows(observationIndices), 0);
        warpingParameters.rows(observationIndices) /= meanParameters(0);
    }
}

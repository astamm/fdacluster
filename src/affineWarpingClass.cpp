#include "affineWarpingClass.h"

unsigned int AffineWarpingFunction::GetNumberOfParameters()
{
    return 2;
}

arma::rowvec AffineWarpingFunction::GetInitialPoint()
{
    return { 1, 0 };
}

arma::mat AffineWarpingFunction::ApplyWarping(const arma::mat &inputGrids,
                                              const arma::mat &warpingParameters)
{
    arma::mat outputGrids(inputGrids.n_rows, inputGrids.n_cols);

    for (unsigned int i = 0;i < inputGrids.n_rows;++i)
        outputGrids.row(i) = inputGrids.row(i) * warpingParameters(i, 0) + warpingParameters(i, 1);

    return outputGrids;
}

void AffineWarpingFunction::SetParameterBounds(const arma::rowvec &warpingOptions,
                                               const arma::mat &inputGrids)
{
    double dl = warpingOptions(0);
    if (dl < 0 || dl > 1)
        Rcpp::stop("The warping option dl for the dilation parameter should be in [0, 1], as the bounds are computed as [1-dl, 1+dl] centered around the unit dilation.");

    double sl = warpingOptions(1);
    double minRange = arma::min(arma::max(inputGrids, 1) - arma::min(inputGrids, 1));

    m_ParameterLowerBounds = { 1 - dl, -sl * minRange};
    m_ParameterUpperBounds = { 1 + dl,  sl * minRange};
}

arma::mat AffineWarpingFunction::GetFinalWarping(const arma::cube &warpingParametersContainer,
                                                 const arma::urowvec &observationMemberships,
                                                 const arma::urowvec &clusterIndices)
{
    unsigned int numberOfObservations = warpingParametersContainer.n_rows;
    unsigned int numberOfParameters = warpingParametersContainer.n_cols;
    unsigned int numberOfIterations = warpingParametersContainer.n_slices;
    arma::mat warpingParameters(numberOfObservations, numberOfParameters, arma::fill::zeros);
    warpingParameters.col(0).ones();
    arma::colvec dilationParameters;
    arma::colvec shiftParameters;

    for (unsigned int i = 0;i < numberOfIterations;++i)
    {
        dilationParameters = warpingParametersContainer.slice(i).col(0);
        shiftParameters = warpingParametersContainer.slice(i).col(1);
        warpingParameters.col(0) %= dilationParameters;
        warpingParameters.col(1) %= dilationParameters;
        warpingParameters.col(1) += shiftParameters;
    }

    this->Normalize(warpingParameters, clusterIndices, observationMemberships);

    return warpingParameters;
}

void AffineWarpingFunction::Normalize(arma::mat &warpingParameters,
                                      const arma::urowvec &clusterIndices,
                                      const arma::urowvec &observationMemberships)
{
    return;

    arma::uvec observationIndices;
    arma::rowvec meanParameters;
    arma::mat clusterValues;

    for (unsigned int i = 0;i < clusterIndices.size();++i)
    {
        observationIndices = arma::find(observationMemberships == clusterIndices(i));
        clusterValues = warpingParameters.rows(observationIndices);

        meanParameters = arma::mean(clusterValues, 0);
        clusterValues.col(0) /= meanParameters(0);
        clusterValues.col(1) -= (meanParameters(1) * clusterValues.col(0));
        warpingParameters.rows(observationIndices) = clusterValues;
    }
}

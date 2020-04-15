#include "affineWarpingClass.h"

unsigned int AffineWarpingFunction::GetNumberOfParameters()
{
    return 2;
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
    unsigned int numberOfParameters = warpingParametersContainer.n_rows;
    unsigned int numberOfObservations = warpingParametersContainer.n_cols;
    unsigned int numberOfIterations = warpingParametersContainer.n_slices;
    arma::mat outputWarpingParameters(numberOfObservations, numberOfParameters, arma::fill::zeros);
    outputWarpingParameters.col(0).ones();
    arma::colvec dilationParameters;
    arma::colvec shiftParameters;

    for (unsigned int i = 0;i < numberOfIterations;++i)
    {
        dilationParameters = warpingParametersContainer.slice(i).col(0);
        shiftParameters = warpingParametersContainer.slice(i).col(1);
        outputWarpingParameters.col(0) %= dilationParameters;
        outputWarpingParameters.col(1) %= dilationParameters;
        outputWarpingParameters.col(1) += shiftParameters;
    }

    arma::uvec observationIndices;
    arma::rowvec meanParameters;
    arma::mat clusterValues;

    for (unsigned int k = 0;k < clusterIndices.size();++k)
    {
        observationIndices = arma::find(observationMemberships == clusterIndices(k));
        meanParameters = arma::mean(outputWarpingParameters.rows(observationIndices), 0);
        clusterValues = outputWarpingParameters.rows(observationIndices);
        clusterValues.col(0) /= meanParameters(0);
        clusterValues.col(1) -= meanParameters(1);
        clusterValues.col(1) /= meanParameters(0);
        outputWarpingParameters.rows(observationIndices) = clusterValues;
    }

    return outputWarpingParameters;
}

void AffineWarpingFunction::Normalize(arma::mat &warpingParameters,
                                      const arma::urowvec &clusterIndices,
                                      const arma::urowvec &observationMemberships)
{
    arma::uvec observationIndices;
    arma::rowvec meanParameters;
    arma::mat clusterValues;

    for (unsigned int i = 0;i < clusterIndices.size();++i)
    {
        observationIndices = arma::find(observationMemberships == clusterIndices(i));
        meanParameters = arma::mean(warpingParameters.rows(observationIndices), 0);

        clusterValues = warpingParameters.rows(observationIndices);
        clusterValues.col(0) /= meanParameters(0);
        clusterValues.col(1) -= meanParameters(1);
        clusterValues.col(1) /= meanParameters(0);
        warpingParameters.rows(observationIndices) = clusterValues;

        // // aggiorno shift e dilation
        // for (unsigned int j = 0;j < observationIndices.size();++j)
        // {
        //     // normalized dilation
        //     warpingParameters(0, observationIndices(j)) = warpingParameters(0, observationIndices(j)) / meanParameters(0);
        //     // normalized shift
        //     warpingParameters(1, observationIndices(j)) = -warpingParameters(0, observationIndices(j)) * meanParameters(1) / meanParameters(0) + warpingParameters(1, observationIndices(j));
        // }
    }
}

double AffineWarpingFunction::GetDissimilarityAfterWarping(const WarpingSet &warpingSet,
                                                           const arma::rowvec &warpingParameters)
{
    return warpingSet.dissimilarityPointer->GetDistance(
            warpingParameters(0) * warpingSet.inputGrid1 + warpingParameters(1),
            warpingSet.inputGrid2,
            warpingSet.inputValues1,
            warpingSet.inputValues2
    );
}

#include "shiftWarpingClass.h"

unsigned int ShiftWarpingFunction::GetNumberOfParameters()
{
    return 1;
}

arma::rowvec ShiftWarpingFunction::GetInitialPoint()
{
    return { 0 };
}

arma::mat ShiftWarpingFunction::ApplyWarping(const arma::mat &inputGrids,
                                             const arma::mat &warpingParameters)
{
    arma::mat outputGrids(inputGrids.n_rows, inputGrids.n_cols);

    for (unsigned int i = 0;i < inputGrids.n_rows;++i)
        outputGrids.row(i) = inputGrids.row(i) + warpingParameters(i, 0);

    return outputGrids;
}

void ShiftWarpingFunction::SetParameterBounds(const arma::rowvec &warpingOptions,
                                              const arma::mat &inputGrids)
{
    double sl = warpingOptions(0);
    double minRange = arma::min(arma::max(inputGrids, 1) - arma::min(inputGrids, 1));

    m_ParameterLowerBounds = { -sl * minRange };
    m_ParameterUpperBounds = {  sl * minRange };
}

arma::mat ShiftWarpingFunction::GetFinalWarping(const arma::cube &warpingParametersContainer,
                                                const arma::urowvec &observationMemberships,
                                                const arma::urowvec &clusterIndices)
{
    unsigned int numberOfObservations = warpingParametersContainer.n_rows;
    unsigned int numberOfParameters = warpingParametersContainer.n_cols;
    unsigned int numberOfIterations = warpingParametersContainer.n_slices;
    arma::mat warpingParameters(numberOfObservations, numberOfParameters, arma::fill::zeros);
    arma::colvec shiftParameters;

    for (unsigned int i = 0;i < numberOfIterations;++i)
    {
        shiftParameters = warpingParametersContainer.slice(i).col(0);
        warpingParameters.col(0) += shiftParameters;
    }

    this->Normalize(warpingParameters, clusterIndices, observationMemberships);

    return warpingParameters;
}

void ShiftWarpingFunction::Normalize(arma::mat &warpingParameters,
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
        warpingParameters.rows(observationIndices) -= meanParameters(0);
    }
}

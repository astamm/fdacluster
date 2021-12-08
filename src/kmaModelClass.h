#ifndef KMAMODELCLASS_H
#define KMAMODELCLASS_H

#include <RcppArmadillo.h>

#include "baseDissimilarityClass.h"
#include "baseWarpingClass.h"
#include "baseCenterClass.h"
#include "baseOptimizerClass.h"

enum ParallelType
{
    ClusterLoop,
    DistanceLoop
};

/// Main class.
/** This class handles loading of the problem and execution of the algorithm.
 */
class KmaModel
{
public:
    KmaModel()
    {
        m_InputGrids.reset();
        m_InputValues.reset();
        m_SeedVector.reset();
        m_WarpingOptions.reset();

        m_NumberOfClusters = 1;
        m_MaximumNumberOfIterations = 100;
        m_NumberOfObservations = 1;
        m_NumberOfDimensions = 1;
        m_NumberOfPoints = 1;
        m_NumberOfThreads = 1;
        m_ParallelMethod = ClusterLoop;

        m_DistanceRelativeTolerance = 1.0e-3;

        m_UseFence = false;
        m_CheckTotalDissimilarity = true;
        m_UseVerbose = true;
        m_ComputeOverallCenter = false;

        std::string m_WarpingMethod = "affine";
        std::string m_CenterMethod = "mean";
        std::string m_DissimilarityMethod = "pearson";
        std::string m_OptimizerMethod = "bobyqa";
    }

    void SetInputData(const arma::mat &grids, const arma::cube &values);
    void SetWarpingMethod(const std::string &val);
    void SetCenterMethod(const std::string &val);
    void SetDissimilarityMethod(const std::string &val);
    void SetOptimizerMethod(const std::string &val);

    void SetSeedVector(const arma::urowvec &val) {m_SeedVector = val;}
    void SetWarpingOptions(const arma::rowvec &val) {m_WarpingOptions = val;}

    void SetNumberOfClusters(const unsigned int &val) {m_NumberOfClusters = val;}
    void SetMaximumNumberOfIterations(const unsigned int &val) {m_MaximumNumberOfIterations = val;}
    void SetNumberOfThreads(const unsigned int &val) {m_NumberOfThreads = val;}
    void SetParallelMethod(const unsigned int &val) {m_ParallelMethod = ParallelType(val);}

    void SetDistanceRelativeTolerance(const double &val) {m_DistanceRelativeTolerance = val;}

    void SetUseFence(const bool &val) {m_UseFence = val;}
    void SetCheckTotalDissimilarity(const bool &val) {m_CheckTotalDissimilarity = val;}
    void SetUseVerbose(const bool &val) {m_UseVerbose = val;}
    void SetComputeOverallCenter(const bool &val) {m_ComputeOverallCenter = val;}

    // Method to get a description of the model.
    void Print(
            const std::string &warpingMethod,
            const std::string &centerMethod,
            const std::string &dissimilarityMethod,
            const std::string &optimizerMethod
    );

    void AlignAndAssignObservations(
            arma::mat &warpingParameters,
            arma::rowvec &observationDistances,
            arma::urowvec &observationMemberships,
            const unsigned int numberOfClusters,
            const arma::urowvec &clusterIndices,
            const arma::mat &warpedGrids,
            const arma::mat &templateGrids,
            const arma::cube &templateValues
    );

    /// Remove warping outliers
    /** It is an optional check that can be activated by the user. After each computation
     *  of best warping parameters, if the computed parameters are flagged as outliers,
     *  they are recomputed with stricter bounds. It is however computationally less
     *  expensive to rerun the kma function with stricter warping options.
     */
    void RunAdaptiveFenceAlgorithm(
            arma::mat &warpingParameters,
            arma::rowvec &observationDistances,
            arma::urowvec &observationMemberships,
            const arma::urowvec &clusterIndices,
            const arma::mat &warpedGrids,
            const arma::mat &templateGrids,
            const arma::cube &templateValues,
            const unsigned int numberOfClusters,
            const double penalizationStep = 0.05
    );

    // Update templates.
    void UpdateTemplates(
            const unsigned int numberOfIterations,
            const arma::urowvec &clusterIndices,
            const arma::urowvec &observationMemberships,
            arma::mat &warpedGrids,
            arma::mat &templateGrids,
            arma::cube &templateValues,
            arma::cube &warpingParametersContainer
    );

    /// Method to execute the algorithm.
    Rcpp::List FitModel();

private:
    arma::mat m_InputGrids;
    arma::cube m_InputValues;
    arma::urowvec m_SeedVector;
    arma::rowvec m_WarpingOptions;

    unsigned int m_NumberOfClusters;
    unsigned int m_MaximumNumberOfIterations;
    unsigned int m_NumberOfObservations;
    unsigned int m_NumberOfDimensions;
    unsigned int m_NumberOfPoints;
    unsigned int m_NumberOfThreads;

    enum ParallelType m_ParallelMethod;

    double m_DistanceRelativeTolerance;

    bool m_UseFence;
    bool m_CheckTotalDissimilarity;
    bool m_UseVerbose;
    bool m_ComputeOverallCenter;

    std::shared_ptr<BaseWarpingFunction> m_WarpingPointer;
    std::shared_ptr<BaseDissimilarityFunction> m_DissimilarityPointer;
    std::shared_ptr<BaseCenterMethod> m_CenterPointer;
    std::shared_ptr<BaseOptimizerFunction> m_OptimizerPointer;
};

#endif /* KMAMODELCLASS_H */

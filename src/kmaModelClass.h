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

        m_NumberOfClusters = 1;
        m_MaximumNumberOfIterations = 100;
        m_NumberOfObservations = 1;
        m_NumberOfDimensions = 1;
        m_NumberOfPoints = 1;
        m_NumberOfThreads = 1;
        m_ParallelMethod = ClusterLoop;

        m_ShiftUpperBound = 0.15;
        m_DilationUpperBound = 0.15;
        m_DistanceRelativeTolerance = 1.0e-3;

        m_UseFence = false;
        m_CheckTotalSimilarity = true;
        m_UseVerbose = true;
        m_ComputeOverallCenter = false;

        std::string m_InterpolationMethod = "linear";
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

    void SetNumberOfClusters(const unsigned int &val) {m_NumberOfClusters = val;}
    void SetMaximumNumberOfIterations(const unsigned int &val) {m_MaximumNumberOfIterations = val;}
    void SetNumberOfThreads(const unsigned int &val) {m_NumberOfThreads = val;}
    void SetParallelMethod(const unsigned int &val) {m_ParallelMethod = ParallelType(val);}

    void SetShiftUpperBound(const double &val) {m_ShiftUpperBound = val;}
    void SetDilationUpperBound(const double &val) {m_DilationUpperBound = val;}
    void SetDistanceRelativeTolerance(const double &val) {m_DistanceRelativeTolerance = val;}

    void SetUseFence(const bool &val) {m_UseFence = val;}
    void SetCheckTotalSimilarity(const bool &val) {m_CheckTotalSimilarity = val;}
    void SetUseVerbose(const bool &val) {m_UseVerbose = val;}
    void SetComputeOverallCenter(const bool &val) {m_ComputeOverallCenter = val;}

    void SetInterpolationMethod(const std::string &val) {m_InterpolationMethod = val;}

    // Method to get a description of the model.
    void Print(
            const std::string &warpingMethod,
            const std::string &centerMethod,
            const std::string &dissimilarityMethod,
            const std::string &optimizerMethod
    );

    // Update templates.
    void UpdateTemplates(
            const arma::mat& x_reg,
            const arma::urowvec& ict,
            const arma::urowvec& labels,
            arma::cube& templates
    );

    /// Method to execute the algorithm.
    Rcpp::List FitModel();

private:
    arma::mat m_InputGrids;
    arma::cube m_InputValues;
    arma::urowvec m_SeedVector;

    unsigned int m_NumberOfClusters;
    unsigned int m_MaximumNumberOfIterations;
    unsigned int m_NumberOfObservations;
    unsigned int m_NumberOfDimensions;
    unsigned int m_NumberOfPoints;
    unsigned int m_NumberOfThreads;

    enum ParallelType m_ParallelMethod;

    double m_ShiftUpperBound;
    double m_DilationUpperBound;
    double m_DistanceRelativeTolerance;

    bool m_UseFence;
    bool m_CheckTotalSimilarity;
    bool m_UseVerbose;
    bool m_ComputeOverallCenter;

    std::string m_InterpolationMethod;

    std::shared_ptr<BaseWarpingFunction> m_WarpingPointer;
    std::shared_ptr<BaseDissimilarityFunction> m_DissimilarityPointer;
    std::shared_ptr<BaseCenterMethod> m_CenterPointer;
    std::shared_ptr<BaseOptimizerFunction> m_OptimizerPointer;
};

#endif /* KMAMODELCLASS_H */

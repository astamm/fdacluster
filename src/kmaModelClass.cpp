#include "kmaModelClass.h"
#include "bobyqaOptimizerClass.h"
#include "noWarpingClass.h"
#include "shiftWarpingClass.h"
#include "dilationWarpingClass.h"
#include "affineWarpingClass.h"
#include "powerWarpingClass.h"
#include "medoidCenterClass.h"
#include "meanCenterClass.h"
#include "pearsonDissimilarityClass.h"
#include "l2DissimilarityClass.h"
#include "powerDissimilarityClass.h"

#include "utilityFunctions.h"
#include "sharedFactoryClass.h"

#include <Rcpp/Benchmark/Timer.h>

#ifdef _OPENMP
#include <omp.h>
#endif

void KmaModel::SetInputData(const arma::mat &grids, const arma::cube &values)
{
  m_InputGrids = grids;
  m_InputValues = values;
  m_NumberOfObservations = values.n_rows;
  m_NumberOfDimensions = values.n_cols;
  m_NumberOfPoints = values.n_slices;

  if (m_InputGrids.n_rows != m_NumberOfObservations)
    Rcpp::stop("The number of observations in the grids does not match the number of observations in the values.");

  if (m_InputGrids.n_cols != m_NumberOfPoints)
    Rcpp::stop("The number of points in the grids does not match the number of points in the values.");
}

void KmaModel::SetWarpingMethod(const std::string &val)
{
  // Warping factory
  SharedFactory<BaseWarpingFunction> warpingFactory;
  warpingFactory.Register<NoWarpingFunction>("none");
  warpingFactory.Register<ShiftWarpingFunction>("shift");
  warpingFactory.Register<DilationWarpingFunction>("dilation");
  warpingFactory.Register<AffineWarpingFunction>("affine");
  warpingFactory.Register<PowerWarpingFunction>("power");

  m_WarpingPointer = warpingFactory.Instantiate(val);

  if (!m_WarpingPointer)
    Rcpp::stop("The warping method is not available.");
}

void KmaModel::SetCenterMethod(const std::string &val)
{
  // Center factory
  SharedFactory<BaseCenterMethod> centerFactory;
  centerFactory.Register<MedoidCenterMethod>("medoid");
  centerFactory.Register<MeanCenterMethod>("mean");

  m_CenterPointer = centerFactory.Instantiate(val);

  if (!m_CenterPointer)
    Rcpp::stop("The center method is not available.");
}

void KmaModel::SetDissimilarityMethod(const std::string &val)
{
  // Dissimilarity factory
  SharedFactory<BaseDissimilarityFunction> dissimilarityFactory;
  dissimilarityFactory.Register<PearsonDissimilarityFunction>("pearson");
  dissimilarityFactory.Register<L2DissimilarityFunction>("l2");
  dissimilarityFactory.Register<PowerDissimilarityFunction>("power");

  m_DissimilarityPointer = dissimilarityFactory.Instantiate(val);

  if (!m_DissimilarityPointer)
    Rcpp::stop("The dissimilarity method is not available.");
}

void KmaModel::SetOptimizerMethod(const std::string &val)
{
  // Optimizer factory
  SharedFactory<BaseOptimizerFunction> optimizerFactory;
  optimizerFactory.Register<BobyqaOptimizerFunction>("bobyqa");

  m_OptimizerPointer = optimizerFactory.Instantiate(val);

  if (!m_OptimizerPointer)
    Rcpp::stop("The optimizer method is not available.");
}

void KmaModel::Print(const std::string &warpingMethod,
                     const std::string &centerMethod,
                     const std::string &dissimilarityMethod,
                     const std::string &optimizerMethod)
{
  Rcpp::Rcout << "Information about the data set:" << std::endl;
  Rcpp::Rcout << " - Number of observations: " << m_NumberOfObservations << std::endl;
  Rcpp::Rcout << " - Number of dimensions: " << m_NumberOfDimensions << std::endl;
  Rcpp::Rcout << " - Number of points: " << m_NumberOfPoints << std::endl;
  Rcpp::Rcout << std::endl;

  Rcpp::Rcout << "Information about cluster initialization:" << std::endl;
  Rcpp::Rcout << " - Number of clusters: " << m_NumberOfClusters << std::endl;
  Rcpp::Rcout << " - Initial seeds for cluster centers: " << m_SeedVector + 1 << std::endl;

  Rcpp::Rcout << "Information about the methods used within the algorithm:" << std::endl;
  Rcpp::Rcout << " - Warping method: " << warpingMethod << std::endl;
  Rcpp::Rcout << " - Center method: " << centerMethod << std::endl;
  Rcpp::Rcout << " - Dissimilarity method: " << dissimilarityMethod << std::endl;
  Rcpp::Rcout << " - Optimization method: " << optimizerMethod << std::endl;
  Rcpp::Rcout << std::endl;

  Rcpp::Rcout << "Information about warping parameter bounds:" << std::endl;
  Rcpp::Rcout << " - Warping options: " << m_WarpingOptions << std::endl;

  Rcpp::Rcout << "Information about convergence criteria:" << std::endl;
  Rcpp::Rcout << " - Maximum number of iterations: " << m_MaximumNumberOfIterations << std::endl;
  Rcpp::Rcout << " - Distance relative tolerance: " << m_DistanceRelativeTolerance << std::endl;
  Rcpp::Rcout << std::endl;

  Rcpp::Rcout << "Information about parallelization setup:" << std::endl;
  Rcpp::Rcout << " - Number of threads: " << m_NumberOfThreads << std::endl;
  Rcpp::Rcout << " - Parallel method: " << m_ParallelMethod << std::endl;
  Rcpp::Rcout << std::endl;

  Rcpp::Rcout << "Other information:" << std::endl;
  Rcpp::Rcout << " - Use fence to robustify: " << m_UseFence << std::endl;
  Rcpp::Rcout << " - Check total dissimilarity: " << m_CheckTotalDissimilarity << std::endl;
  Rcpp::Rcout << " - Compute overall center: " << m_ComputeOverallCenter << std::endl;
  Rcpp::Rcout << std::endl;
}

void KmaModel::AlignAndAssignObservation(arma::mat &warpingParameters,
                                         arma::rowvec &observationDistances,
                                         arma::urowvec &observationMemberships,
                                         const unsigned int observationIndex,
                                         const unsigned int numberOfClusters,
                                         const arma::urowvec &clusterIndices,
                                         const arma::mat &warpedGrids,
                                         const arma::mat &templateGrids,
                                         const arma::cube &templateValues)
{
  unsigned int numberOfParameters = warpingParameters.n_cols;
  arma::rowvec workingObservationDistances(numberOfClusters);
  arma::mat workingParameterValues(numberOfClusters, numberOfParameters);
  arma::rowvec startingParameters(numberOfParameters);
  WarpingSet warpingSet;

  // Compute warping parameters for each template
  for (unsigned int i = 0;i < numberOfClusters;++i)
  {
    warpingSet = m_WarpingPointer->SetInputData(
      warpedGrids.row(observationIndex),
      templateGrids.row(i),
      m_InputValues.row(observationIndex),
      templateValues.row(i),
      m_DissimilarityPointer
    );

    workingObservationDistances(i) = m_OptimizerPointer->Optimize(
      startingParameters,
      m_WarpingPointer,
      warpingSet
    );

    workingParameterValues.row(i) = startingParameters;
  }

  observationDistances(observationIndex) = workingObservationDistances.min();
  unsigned int assignedTemplateIndex = arma::index_min(workingObservationDistances);
  observationMemberships(observationIndex) = clusterIndices(assignedTemplateIndex);
  warpingParameters.row(observationIndex) = workingParameterValues.row(assignedTemplateIndex);
}

void KmaModel::RunAdaptiveFenceAlgorithm(arma::mat &warpingParameters,
                                         arma::rowvec &observationDistances,
                                         arma::urowvec &observationMemberships,
                                         const arma::urowvec &clusterIndices,
                                         const arma::mat &warpedGrids,
                                         const arma::mat &templateGrids,
                                         const arma::cube &templateValues,
                                         const unsigned int numberOfClusters,
                                         const unsigned int maximumNumberOfIterations)
{
  unsigned int numberOfParameters = warpingParameters.n_cols;
  arma::vec quantileOrders = { 0.25, 0.75 };
  arma::mat quantileValues;
  arma::mat reasonableBounds;
  arma::urowvec outlierIndices;
  arma::urowvec workingIndices;
  unsigned int runningIteration = 0;
  bool continueLoop = true;

  while (continueLoop)
  {
    quantileValues = arma::quantile(warpingParameters, quantileOrders);
    reasonableBounds = quantileValues;
    reasonableBounds.row(0) -= 1.5 * (quantileValues.row(1) - quantileValues.row(0));
    reasonableBounds.row(1) += 1.5 * (quantileValues.row(1) - quantileValues.row(0));

    outlierIndices.reset();

    for (unsigned int i = 0;i < numberOfParameters;++i)
    {
      workingIndices = arma::find(warpingParameters.col(i) < reasonableBounds(0, i) || warpingParameters.col(i) > reasonableBounds(1, i));
      outlierIndices = arma::join_horiz(outlierIndices, workingIndices);
    }

    outlierIndices = arma::unique(outlierIndices);

    if (outlierIndices.size() == 0)
    {
      continueLoop = false;
      continue;
    }

    // Redo optimization foe each observation with outliers
    m_WarpingPointer->SetParameterBounds(reasonableBounds);

#ifdef _OPENMP
#pragma omp parallel for num_threads(m_NumberOfThreads)
#endif

    for (unsigned int i = 0;i < outlierIndices.size();++i)
    {
      unsigned int observationIndex = outlierIndices(i);

      this->AlignAndAssignObservation(
          warpingParameters,
          observationDistances,
          observationMemberships,
          observationIndex,
          numberOfClusters,
          clusterIndices,
          warpedGrids,
          templateGrids,
          templateValues
      );
    }

    ++runningIteration;

    if (runningIteration >= maximumNumberOfIterations)
      continueLoop = false;
  }

  return;
}

void KmaModel::UpdateTemplates(const arma::mat& warpedGrids,
                               const arma::urowvec& clusterIndices,
                               const arma::urowvec& observationMemberships,
                               arma::mat& templateGrids,
                               arma::cube& templateValues)
{
  // switch to choose how to parallelize
  // case ClusterLoop: each thread one cluster
  // case DistanceLoop: each cluster all the threads (available only with medoid)

  switch(m_ParallelMethod)
  {
  case ClusterLoop:

#ifdef _OPENMP
#pragma omp parallel for num_threads(m_NumberOfThreads)
#endif

    for (unsigned int i = 0;i < clusterIndices.size();++i)
    {
      arma::uvec selectedObservations = arma::find(observationMemberships == clusterIndices(i));

      CenterType centerComputer = m_CenterPointer->GetCenter(
        warpedGrids.rows(selectedObservations),
        GetObservations(m_InputValues, selectedObservations),
        m_DissimilarityPointer
      );

      templateGrids.row(i) = centerComputer.centerGrid;
      templateValues.row(i) = centerComputer.centerValues;
    }

    break;

  case DistanceLoop:

    arma::uvec selectedObservations;
    CenterType centerComputer;

    for (unsigned int i = 0;i < clusterIndices.size();++i)
    {
      selectedObservations = arma::find(observationMemberships == clusterIndices(i));

      centerComputer = m_CenterPointer->GetCenter(
        warpedGrids.rows(selectedObservations),
        GetObservations(m_InputValues, selectedObservations),
        m_DissimilarityPointer,
        m_NumberOfThreads
      );

      templateGrids.row(i) = centerComputer.centerGrid;
      templateValues.row(i) = centerComputer.centerValues;
    }

    break;
  }
}

Rcpp::List KmaModel::FitModel()
{
  Rcpp::Timer timer;
  timer.step("start execution");

  //
  // initial templates
  //
  arma::mat templateGrids(m_NumberOfClusters, m_NumberOfPoints);
  arma::cube templateValues(m_NumberOfClusters, m_NumberOfDimensions, m_NumberOfPoints);

  for (unsigned int i = 0;i < m_NumberOfClusters;++i)
  {
    templateGrids.row(i) = m_InputGrids.row(m_SeedVector(i));
    templateValues.row(i) = m_InputValues.row(m_SeedVector(i));
  }

  // Initialize containers for storing
  // template grids and values at each iteration
  arma::cube templateGridsContainer(m_NumberOfClusters, m_NumberOfPoints, m_MaximumNumberOfIterations);
  arma::field<arma::cube> templateValuesContainer(m_MaximumNumberOfIterations);

  templateGridsContainer.slice(0) = templateGrids;
  templateValuesContainer(0) = templateValues;

  //
  // compute center_origin (to be fixed with new centers)
  //

  CenterType overallCenter;

  if (m_ComputeOverallCenter)
    overallCenter = m_CenterPointer->GetCenter(
      m_InputGrids,
      m_InputValues,
      m_DissimilarityPointer
    );

  //
  // WHILE equipment
  //
  unsigned int numberOfParameters = m_WarpingPointer->GetNumberOfParameters();
  unsigned int numberOfIterations = 0;
  unsigned int numberOfClusters = m_NumberOfClusters;
  arma::rowvec observationDistances(m_NumberOfObservations, arma::fill::ones);
  arma::rowvec oldObservationDistances(m_NumberOfObservations, arma::fill::zeros);
  arma::urowvec observationMemberships(m_NumberOfObservations, arma::fill::ones);
  arma::urowvec oldObservationMemberships(m_NumberOfObservations, arma::fill::zeros);
  arma::urowvec clusterIndices = arma::linspace<arma::urowvec>(0, m_NumberOfClusters - 1, m_NumberOfClusters);
  arma::mat warpedGrids = m_InputGrids;
  arma::mat warpingParameters(m_NumberOfObservations, numberOfParameters);
  arma::cube warpingParametersContainer(m_NumberOfObservations, numberOfParameters, m_MaximumNumberOfIterations);
  bool distanceCondition = true;
  bool membershipCondition = true;
  bool iterationCondition = true;
  bool totalDissimilarityCondition = true;

  timer.step("seeds and original center");

  if (m_UseVerbose)
    Rcpp::Rcout << "Running k-centroid algorithm:" << std::endl;

  while (distanceCondition && membershipCondition && iterationCondition && totalDissimilarityCondition)
  {
    ++numberOfIterations;
    iterationCondition = numberOfIterations < m_MaximumNumberOfIterations;

    if (m_UseVerbose)
      Rcpp::Rcout << " * Iteration #" << numberOfIterations << std::endl;

    oldObservationDistances = observationDistances;
    oldObservationMemberships = observationMemberships;

    m_WarpingPointer->SetParameterBounds(m_WarpingOptions, warpedGrids);

#ifdef _OPENMP
#pragma omp parallel for num_threads(m_NumberOfThreads)
#endif

    for (unsigned int i = 0;i < m_NumberOfObservations;++i)
    {
      this->AlignAndAssignObservation(
          warpingParameters,
          observationDistances,
          observationMemberships,
          i,
          numberOfClusters,
          clusterIndices,
          warpedGrids,
          templateGrids,
          templateValues
      );
    }

    distanceCondition = arma::any(arma::abs(observationDistances - oldObservationDistances) > m_DistanceRelativeTolerance * observationDistances);
    membershipCondition = arma::any(observationMemberships != oldObservationMemberships) || (m_NumberOfClusters == 1);

    timer.step( "warping "+ std::to_string(numberOfIterations) );

    // Update current template list
    clusterIndices = arma::unique(observationMemberships);
    numberOfClusters = clusterIndices.size();

    // Print cluster sizes
    if (m_UseVerbose)
    {
      std::map<unsigned int,unsigned int> labelCounts = tableCpp(observationMemberships);
      for (auto it = labelCounts.cbegin();it != labelCounts.cend();++it)
        Rcpp::Rcout <<"   - Size of cluster #" << it->first << ": " << it->second << std::endl;
    }

    if (m_UseFence)
      this->RunAdaptiveFenceAlgorithm(
          warpingParameters,
          observationDistances,
          observationMemberships,
          clusterIndices,
          warpedGrids,
          templateGrids,
          templateValues,
          numberOfClusters
      );

    // Normalization
    m_WarpingPointer->Normalize(warpingParameters, clusterIndices, observationMemberships);

    // Store parameter values in container
    warpingParametersContainer.slice(numberOfIterations - 1) = warpingParameters;

    // Update individual warped grids
    warpedGrids = m_WarpingPointer->ApplyWarping(warpedGrids, warpingParameters);

    timer.step( "fence/norm/update "+ std::to_string(numberOfIterations) );

    // Compute new templates
    templateGrids.set_size(numberOfClusters, m_NumberOfPoints);
    templateValues.set_size(numberOfClusters, m_NumberOfDimensions, m_NumberOfPoints);

    this->UpdateTemplates(
        warpedGrids,
        clusterIndices,
        observationMemberships,
        templateGrids,
        templateValues
    );

    templateGridsContainer.slice(numberOfIterations) = templateGrids;
    templateValuesContainer(numberOfIterations) = templateValues;

    // Check total dissimilarity
    if (m_CheckTotalDissimilarity)
    {
      double totalDissimilarity = arma::sum(observationDistances);
      double oldTotalDissimilarity = arma::sum(oldObservationDistances);

      // if total distance increased
      if (oldTotalDissimilarity <= totalDissimilarity)
      {
        totalDissimilarityCondition = false;
        templateGrids = templateGridsContainer.slice(numberOfIterations - 1);
        templateValues = templateValuesContainer(numberOfIterations - 1);
        observationDistances = oldObservationDistances;
        observationMemberships = oldObservationMemberships;
        --numberOfIterations;
      }
    }

    timer.step( "newtemplates "+ std::to_string(numberOfIterations) );
  }

  warpingParametersContainer.resize(m_NumberOfObservations, numberOfParameters, numberOfIterations);

  // Compute final warping
  arma::mat finalWarpingParameters = m_WarpingPointer->GetFinalWarping(
    warpingParametersContainer,
    observationMemberships,
    clusterIndices
  );

  // Convert cube to field for conversion to List in R
  arma::field<arma::mat> listOfEstimatedParameters(numberOfIterations);
  arma::field<arma::mat> listOfTemplateGrids(numberOfIterations + 1);

  listOfTemplateGrids(0) = templateGridsContainer.slice(0);
  for (unsigned int k = 0;k < numberOfIterations;++k)
  {
    listOfEstimatedParameters(k) = warpingParametersContainer.slice(k);
    listOfTemplateGrids(k + 1) = templateGridsContainer.slice(k + 1);
  }

  templateValuesContainer = templateValuesContainer.rows(0, numberOfIterations);

  timer.step( "output ");

  // Convert arma vector types to Rcpp::NumericVector
  // since arma types consistently converts to matrix in R
  Rcpp::NumericVector outputSeedVector = FormatVectorForOutput(m_SeedVector);
  Rcpp::NumericVector outputDistancesToOverallCenter = FormatVectorForOutput(overallCenter.distancesToCenter);
  Rcpp::NumericVector outputObservationMemberships = FormatVectorForOutput(observationMemberships + 1);
  Rcpp::NumericVector outputObservationDistances = FormatVectorForOutput(observationDistances);

  return Rcpp::List::create(
    Rcpp::Named("x")                           = m_InputGrids,
    Rcpp::Named("y")                           = m_InputValues,
    Rcpp::Named("seeds")                       = outputSeedVector,
    Rcpp::Named("iterations")                  = numberOfIterations,
    Rcpp::Named("n_clust")                     = m_NumberOfClusters,
    Rcpp::Named("overall_center_grid")         = overallCenter.centerGrid,
    Rcpp::Named("overall_center_values")       = overallCenter.centerValues,
    Rcpp::Named("distances_to_overall_center") = outputDistancesToOverallCenter,
    Rcpp::Named("x_final")                     = warpedGrids,
    Rcpp::Named("n_clust_final")               = clusterIndices.size(),
    Rcpp::Named("x_centers_final")             = templateGrids,
    Rcpp::Named("y_centers_final")             = templateValues,
    Rcpp::Named("template_grids")              = listOfTemplateGrids,
    Rcpp::Named("template_values")             = templateValuesContainer,
    Rcpp::Named("labels")                      = outputObservationMemberships,
    Rcpp::Named("final_dissimilarity")         = outputObservationDistances,
    Rcpp::Named("parameters_list")             = listOfEstimatedParameters,
    Rcpp::Named("parameters")                  = finalWarpingParameters,
    Rcpp::Named("timer")                       = timer
  );
}

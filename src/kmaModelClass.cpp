#include "kmaModelClass.h"
#include "bobyqaOptimizerClass.h"
#include "noWarpingClass.h"
#include "shiftWarpingClass.h"
#include "dilationWarpingClass.h"
#include "affineWarpingClass.h"
#include "medoidCenterClass.h"
#include "meanCenterClass.h"
#include "pearsonDissimilarityClass.h"
#include "l2DissimilarityClass.h"

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

void KmaModel::AlignAndAssignObservations(arma::mat &warpingParameters,
                                          arma::rowvec &observationDistances,
                                          arma::urowvec &observationMemberships,
                                          const unsigned int numberOfClusters,
                                          const arma::urowvec &clusterIndices,
                                          const arma::mat &warpedGrids,
                                          const arma::mat &templateGrids,
                                          const arma::cube &templateValues)
{
  unsigned int numberOfParameters = warpingParameters.n_cols;

#ifdef _OPENMP
#pragma omp parallel for num_threads(m_NumberOfThreads)
#endif

  for (unsigned int j = 0;j < m_NumberOfObservations;++j)
  {
    arma::rowvec workingObservationDistances(numberOfClusters);
    arma::mat workingParameterValues(numberOfClusters, numberOfParameters);
    arma::rowvec startingParameters(numberOfParameters);
    WarpingSet warpingSet;

    // Compute warping parameters for each template
    for (unsigned int i = 0;i < numberOfClusters;++i)
    {
      warpingSet = m_WarpingPointer->SetInputData(
        warpedGrids.row(j),
        templateGrids.row(i),
        m_InputValues.row(j),
        templateValues.row(i),
        m_DissimilarityPointer
      );

      workingObservationDistances(i) = m_OptimizerPointer->AlignToTemplate(
        startingParameters,
        warpedGrids.row(j),
        m_InputValues.row(j),
        templateGrids.row(i),
        templateValues.row(i),
        m_DissimilarityPointer,
        m_WarpingPointer
      );

      workingParameterValues.row(i) = startingParameters;
    }

    observationDistances(j) = workingObservationDistances.min();
    unsigned int assignedTemplateIndex = arma::index_min(workingObservationDistances);
    observationMemberships(j) = clusterIndices(assignedTemplateIndex);
    warpingParameters.row(j) = workingParameterValues.row(assignedTemplateIndex);
  }
}

void KmaModel::RunAdaptiveFenceAlgorithm(arma::mat &warpingParameters,
                                         arma::rowvec &observationDistances,
                                         arma::urowvec &observationMemberships,
                                         const arma::urowvec &clusterIndices,
                                         const arma::mat &warpedGrids,
                                         const arma::mat &templateGrids,
                                         const arma::cube &templateValues,
                                         const unsigned int numberOfClusters,
                                         const double penalizationStep)
{
  unsigned int numberOfParameters = warpingParameters.n_cols;
  arma::vec quantileOrders = { 0.25, 0.75 };
  arma::mat quantileValues;
  arma::mat fenceValues;
  bool continueLoop = true;
  double penalizationWeight = 0.0;

  while (continueLoop)
  {
    quantileValues = arma::quantile(warpingParameters, quantileOrders);

    fenceValues = quantileValues;
    fenceValues.row(0) -= 1.5 * (quantileValues.row(1) - quantileValues.row(0));
    fenceValues.row(1) += 1.5 * (quantileValues.row(1) - quantileValues.row(0));

    continueLoop = false;

    for (unsigned int i = 0;i < numberOfParameters;++i)
    {
      if (arma::any(warpingParameters.col(i) < fenceValues(0, i)))
      {
        continueLoop = true;
        break;
      }

      if (arma::any(warpingParameters.col(i) > fenceValues(1, i)))
      {
        continueLoop = true;
        break;
      }
    }

    if (!continueLoop)
      continue;

    penalizationWeight += penalizationStep;

    if (penalizationWeight > 1.0)
    {
      continueLoop = false;
      continue;
    }

    m_OptimizerPointer->SetPenalizationWeight(penalizationWeight);

    this->AlignAndAssignObservations(
        warpingParameters,
        observationDistances,
        observationMemberships,
        numberOfClusters,
        clusterIndices,
        warpedGrids,
        templateGrids,
        templateValues
    );
  }
}

void KmaModel::UpdateTemplates(const unsigned int numberOfIterations,
                               const arma::urowvec &clusterIndices,
                               const arma::urowvec &observationMemberships,
                               arma::mat &warpedGrids,
                               arma::mat &templateGrids,
                               arma::cube &templateValues,
                               arma::cube &warpingParametersContainer)
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

      arma::rowvec warpingParameters;

      m_OptimizerPointer->CenterTemplate(
          warpingParameters,
          centerComputer.centerGrid,
          centerComputer.centerValues,
          m_InputGrids.rows(selectedObservations),
          GetObservations(m_InputValues, selectedObservations),
          m_DissimilarityPointer,
          m_WarpingPointer
      );

      templateGrids.row(i) = m_WarpingPointer->ApplyWarping(centerComputer.centerGrid, warpingParameters);
      templateValues.row(i) = centerComputer.centerValues;

      for (unsigned int j = 0;j < selectedObservations.size();++j)
      {
        unsigned int observationIndex = selectedObservations(j);
        warpedGrids.row(observationIndex) = m_WarpingPointer->ApplyWarping(warpedGrids.row(observationIndex), warpingParameters);
        warpingParametersContainer.slice(2 * (numberOfIterations - 1) + 1).row(observationIndex) = warpingParameters;
      }
    }

    break;

  case DistanceLoop:

    arma::uvec selectedObservations;
    CenterType centerComputer;
    arma::rowvec warpingParameters;

    for (unsigned int i = 0;i < clusterIndices.size();++i)
    {
      selectedObservations = arma::find(observationMemberships == clusterIndices(i));

      centerComputer = m_CenterPointer->GetCenter(
        warpedGrids.rows(selectedObservations),
        GetObservations(m_InputValues, selectedObservations),
        m_DissimilarityPointer,
        m_NumberOfThreads
      );

      m_OptimizerPointer->CenterTemplate(
          warpingParameters,
          centerComputer.centerGrid,
          centerComputer.centerValues,
          m_InputGrids.rows(selectedObservations),
          GetObservations(m_InputValues, selectedObservations),
          m_DissimilarityPointer,
          m_WarpingPointer
      );

      templateGrids.row(i) = m_WarpingPointer->ApplyWarping(centerComputer.centerGrid, warpingParameters);
      templateValues.row(i) = centerComputer.centerValues;

      for (unsigned int j = 0;j < selectedObservations.size();++j)
      {
        unsigned int observationIndex = selectedObservations(j);
        warpedGrids.row(observationIndex) = m_WarpingPointer->ApplyWarping(warpedGrids.row(observationIndex), warpingParameters);
        warpingParametersContainer.slice(2 * (numberOfIterations - 1) + 1).row(observationIndex) = warpingParameters;
      }
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
  arma::field<arma::mat> templateGridsContainer(m_MaximumNumberOfIterations);
  arma::field<arma::cube> templateValuesContainer(m_MaximumNumberOfIterations);

  templateGridsContainer(0) = templateGrids;
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
  arma::cube warpingParametersContainer(m_NumberOfObservations, numberOfParameters, 2 * m_MaximumNumberOfIterations);
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
      Rcpp::Rcout << " - Iteration #" << numberOfIterations << std::endl;

    oldObservationDistances = observationDistances;
    oldObservationMemberships = observationMemberships;

    m_WarpingPointer->SetParameterBounds(m_WarpingOptions, warpedGrids);
    m_OptimizerPointer->SetPenalizationWeight(0.0);

    this->AlignAndAssignObservations(
        warpingParameters,
        observationDistances,
        observationMemberships,
        numberOfClusters,
        clusterIndices,
        warpedGrids,
        templateGrids,
        templateValues
    );

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
        Rcpp::Rcout <<"   * Size of cluster #" << it->first << ": " << it->second << std::endl;
    }

    // Normalization
    m_WarpingPointer->Normalize(warpingParameters, clusterIndices, observationMemberships);

    // Store parameter values in container
    warpingParametersContainer.slice(2 * (numberOfIterations - 1)) = warpingParameters;

    // Update individual warped grids
    warpedGrids = m_WarpingPointer->ApplyWarping(warpedGrids, warpingParameters);

    timer.step( "fence/norm/update "+ std::to_string(numberOfIterations) );

    // Compute new templates
    templateGrids.set_size(numberOfClusters, m_NumberOfPoints);
    templateValues.set_size(numberOfClusters, m_NumberOfDimensions, m_NumberOfPoints);

    this->UpdateTemplates(
        numberOfIterations,
        clusterIndices,
        observationMemberships,
        warpedGrids,
        templateGrids,
        templateValues,
        warpingParametersContainer
    );

    templateGridsContainer(numberOfIterations) = templateGrids;
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
        templateGrids = templateGridsContainer(numberOfIterations - 1);
        templateValues = templateValuesContainer(numberOfIterations - 1);
        observationDistances = oldObservationDistances;
        observationMemberships = oldObservationMemberships;
        --numberOfIterations;
      }
    }

    timer.step( "newtemplates "+ std::to_string(numberOfIterations) );
  }

  if (m_UseVerbose)
  {
    Rcpp::Rcout << std::endl;
    Rcpp::Rcout << "Active stopping criteria:" << std::endl;

    if (!distanceCondition)
      Rcpp::Rcout << " - Individual distances did not change within the user-defined relative tolerance (" << m_DistanceRelativeTolerance << ")." << std::endl;

    if (!membershipCondition)
      Rcpp::Rcout << " - Memberships did not change." << std::endl;

    if (!iterationCondition)
      Rcpp::Rcout << " - Maximum number of iterations reached (" << m_MaximumNumberOfIterations << ")." << std::endl;

    if (m_CheckTotalDissimilarity && !totalDissimilarityCondition)
      Rcpp::Rcout << " - The total dissimilarity did not decrease." << std::endl;
  }

  warpingParametersContainer.resize(m_NumberOfObservations, numberOfParameters, 2 * numberOfIterations);

  // Compute final warping
  arma::mat finalWarpingParameters = m_WarpingPointer->GetFinalWarping(
    warpingParametersContainer,
    observationMemberships,
    clusterIndices
  );

  timer.step( "output ");

  // Convert arma vector types to Rcpp::NumericVector
  // since arma types consistently converts to matrix in R
  Rcpp::NumericVector outputSeedVector = FormatVectorForOutput(m_SeedVector);
  Rcpp::NumericVector outputDistancesToOverallCenter = FormatVectorForOutput(overallCenter.distancesToCenter);
  Rcpp::NumericVector outputObservationMemberships = FormatVectorForOutput(observationMemberships + 1);
  Rcpp::NumericVector outputObservationDistances = FormatVectorForOutput(observationDistances);

  // Convert containers into Rcpp::List for conversion to list in R
  Rcpp::List listOfEstimatedParameters(numberOfIterations);
  Rcpp::List listOfTemplateGrids(numberOfIterations + 1);
  Rcpp::List listOfTemplateValues(numberOfIterations + 1);

  listOfTemplateGrids[0]  = templateGridsContainer(0);
  listOfTemplateValues[0] = templateValuesContainer(0);

  for (unsigned int k = 0;k < numberOfIterations;++k)
  {
    listOfEstimatedParameters[k] = warpingParametersContainer.slice(k);
    listOfTemplateGrids[k + 1]   = templateGridsContainer(k + 1);
    listOfTemplateValues[k + 1]  = templateValuesContainer(k + 1);
  }

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
    Rcpp::Named("template_values")             = listOfTemplateValues,
    Rcpp::Named("labels")                      = outputObservationMemberships,
    Rcpp::Named("final_dissimilarity")         = outputObservationDistances,
    Rcpp::Named("parameters_list")             = listOfEstimatedParameters,
    Rcpp::Named("parameters")                  = finalWarpingParameters,
    Rcpp::Named("timer")                       = timer
  );
}

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
#include "fenceAlgorithm.h"

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

  Rcpp::Rcout << "Information about cluster initialization:" << std::endl;
  Rcpp::Rcout << " - Number of clusters: " << m_NumberOfClusters << std::endl;
  Rcpp::Rcout << " - Initial seeds for cluster centers: " << m_SeedVector + 1 << std::endl;

  Rcpp::Rcout << "Information about the methods used within the algorithm:" << std::endl;
  Rcpp::Rcout << " - Warping method: " << warpingMethod << std::endl;
  Rcpp::Rcout << " - Center method: " << centerMethod << std::endl;
  Rcpp::Rcout << " - Dissimilarity method: " << dissimilarityMethod << std::endl;
  Rcpp::Rcout << " - Optimization method: " << optimizerMethod << std::endl;

  Rcpp::Rcout << "Information about warping parameter bounds:" << std::endl;
  Rcpp::Rcout << " - Warping options: " << m_WarpingOptions << std::endl;

  Rcpp::Rcout << "Information about convergence criteria:" << std::endl;
  Rcpp::Rcout << " - Maximum number of iterations: " << m_MaximumNumberOfIterations << std::endl;
  Rcpp::Rcout << " - Distance relative tolerance: " << m_DistanceRelativeTolerance << std::endl;

  Rcpp::Rcout << "Information about parallelization setup:" << std::endl;
  Rcpp::Rcout << " - Number of threads: " << m_NumberOfThreads << std::endl;
  Rcpp::Rcout << " - Parallel method: " << m_ParallelMethod << std::endl;

  Rcpp::Rcout << "Other information:" << std::endl;
  Rcpp::Rcout << " - Use fence to robustify: " << m_UseFence << std::endl;
  Rcpp::Rcout << " - Check total similarity: " << m_CheckTotalSimilarity << std::endl;
  Rcpp::Rcout << " - Compute overall center: " << m_ComputeOverallCenter << std::endl;
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

  arma::uvec selectedObservations;
  CenterType centerComputer;

  switch(m_ParallelMethod)
  {
  case ClusterLoop:

#ifdef _OPENMP
#pragma omp parallel for num_threads(m_NumberOfThreads)
#endif

    for (unsigned int i = 0;i < clusterIndices.size();++i)
    {
      selectedObservations = arma::find(observationMemberships == clusterIndices(i));

      centerComputer = m_CenterPointer->GetCenter(
        warpedGrids.rows(selectedObservations),
        GetObservations(m_InputValues, selectedObservations),
        m_DissimilarityPointer
      );

      templateGrids.row(i) = centerComputer.centerGrid;
      templateValues.tube(arma::span(i), arma::span::all) = centerComputer.centerValues;
    }

    break;

  case DistanceLoop:

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
      templateValues.tube(arma::span(i), arma::span::all) = centerComputer.centerValues;

      if (m_UseVerbose)
        Rcpp::Rcout << "Template num. " << i << " updated." << std::endl;
    }

    break;
  }
}

Rcpp::List KmaModel::FitModel()
{
  if (m_UseVerbose)
    Rcpp::Rcout << "Start execution." << std::endl;

  Rcpp::Timer timer;
  timer.step("start execution");

  //
  // initial templates
  //
  if (m_UseVerbose)
    Rcpp::Rcout << "Compute initial templates: ";

  arma::mat templateGrids(m_NumberOfClusters, m_NumberOfPoints);
  arma::cube templateValues(m_NumberOfClusters, m_NumberOfDimensions, m_NumberOfPoints);

  for (unsigned int i = 0;i < m_NumberOfClusters;++i)
  {
    templateGrids.row(i) = m_InputGrids.row(m_SeedVector(i));
    templateValues.tube(arma::span(i), arma::span::all) = GetObservation(m_InputValues, m_SeedVector(i));
  }

  // Initialize containers for storing
  // template grids and values at each iteration
  arma::cube templateGridsContainer(m_NumberOfClusters, m_NumberOfPoints, m_MaximumNumberOfIterations);
  arma::field<arma::cube> templateValuesContainer(m_MaximumNumberOfIterations);

  templateGridsContainer.slice(0) = templateGrids;
  templateValuesContainer(0) = templateValues;

  if (m_UseVerbose)
    Rcpp::Rcout << "Done." << std::endl;

  //
  // compute center_origin (to be fixed with new centers)
  //

  CenterType overallCenter;

  if (m_ComputeOverallCenter)
  {
    if (m_UseVerbose)
      Rcpp::Rcout << "Compute overall center: ";

    overallCenter = m_CenterPointer->GetCenter(m_InputGrids, m_InputValues, m_DissimilarityPointer);

    if (m_UseVerbose)
      Rcpp::Rcout << "Done" << std::endl;
  }

  //
  // WHILE equipment
  //
  if (m_UseVerbose)
    Rcpp::Rcout << "Start while iteration" << std::endl;

  // distances to template at each iteration
  arma::rowvec observationDistances(m_NumberOfObservations, arma::fill::ones);
  arma::rowvec oldObservationDistances(m_NumberOfObservations, arma::fill::zeros);

  // warping parameter container for storing parameter estimates at each iteration
  unsigned int numberOfParameters = m_WarpingPointer->GetNumberOfParameters();
  arma::cube warpingParametersContainer(m_NumberOfObservations, numberOfParameters, m_MaximumNumberOfIterations);
  arma::mat warpedGrids = m_InputGrids;

  // observation memberships at each iteration
  arma::urowvec observationMemberships(m_NumberOfObservations, arma::fill::ones);
  arma::urowvec oldObservationMemberships(m_NumberOfObservations, arma::fill::zeros);

  // cluster indices
  arma::urowvec clusterIndices = arma::linspace<arma::urowvec>(0, m_NumberOfClusters - 1, m_NumberOfClusters);
  unsigned int iter = 0;

  timer.step("seeds and original center");

  bool distanceCondition = arma::any(arma::abs(observationDistances - oldObservationDistances) > m_DistanceRelativeTolerance * observationDistances);
  bool membershipCondition = arma::any(observationMemberships != oldObservationMemberships) || (m_NumberOfClusters == 1);
  bool iterationCondition = iter < m_MaximumNumberOfIterations;
  bool totalDissimilarityCondition = true;

  arma::mat warpingParameters(m_NumberOfObservations, numberOfParameters);
  unsigned int numberOfClusters = m_NumberOfClusters;

  while (distanceCondition && membershipCondition && iterationCondition && totalDissimilarityCondition)
  {
    ++iter;
    iterationCondition = iter < m_MaximumNumberOfIterations;

    if (m_UseVerbose)
      Rcpp::Rcout << "Iteration #" << iter << std::endl;

    oldObservationDistances = observationDistances;
    oldObservationMemberships = observationMemberships;

    if (m_UseVerbose)
      Rcpp::Rcout << "Set bound: ";

    m_WarpingPointer->SetParameterBounds(m_WarpingOptions, warpedGrids);

    if (m_UseVerbose)
      Rcpp::Rcout << "Done" << std::endl;

    // compute best warping parameters and assign new labels
    if (m_UseVerbose)
      Rcpp::Rcout << iter << ". Compute best warping: " << std::endl;

    arma::rowvec workingObservationDistances(numberOfClusters);
    arma::mat workingParameterValues(numberOfClusters, numberOfParameters);
    arma::rowvec startingParameters(numberOfParameters);
    arma::rowvec workingWarpedGrid;
    arma::rowvec workingTemplateGrid;
    arma::mat workingValues;
    arma::mat workingTemplateValues;
    WarpingSet warpingSet;

#ifdef _OPENMP
#pragma omp parallel for num_threads(m_NumberOfThreads)
#endif

    for (unsigned int i = 0;i < m_NumberOfObservations;++i)
    {
      workingWarpedGrid = warpedGrids.row(i);
      workingValues = GetObservation(m_InputValues, i);

      // Compute warping parameters for each template
      for (unsigned int j = 0;j < numberOfClusters;++j)
      {
        workingTemplateGrid = templateGrids.row(j);
        workingTemplateValues = templateValues.tube(arma::span(j), arma::span::all);

        warpingSet = m_WarpingPointer->SetInputData(
          workingWarpedGrid,
          workingTemplateGrid,
          workingValues,
          workingTemplateValues,
          m_DissimilarityPointer
        );

        auto costFunction = [this, &warpingSet] (const arma::rowvec &arg)
        {
          return this->m_WarpingPointer->GetDissimilarityAfterWarping(warpingSet, arg);
        };

        workingObservationDistances(j) = m_OptimizerPointer->Optimize(
          startingParameters,
          m_WarpingPointer,
          costFunction
        );

        workingParameterValues.row(j) = startingParameters;
      }

      observationDistances(i) = workingObservationDistances.min();
      unsigned int assignedTemplateIndex = arma::index_min(workingObservationDistances);
      observationMemberships(i) = clusterIndices(assignedTemplateIndex);
      warpingParameters.row(i) = workingParameterValues.row(assignedTemplateIndex);
    }

    if (m_UseVerbose)
      Rcpp::Rcout << "Done" << std::endl;

    distanceCondition = arma::any(arma::abs(observationDistances - oldObservationDistances) > m_DistanceRelativeTolerance * observationDistances);
    membershipCondition = arma::any(observationMemberships != oldObservationMemberships) || (m_NumberOfClusters == 1);

    timer.step( "warping "+ std::to_string(iter) );

    // Update current template list
    clusterIndices = arma::unique(observationMemberships);
    numberOfClusters = clusterIndices.size();

    //PRINT
    if (m_UseVerbose)
    {
      clusterIndices.print(Rcpp::Rcout, "Cluster indices: ");
      std::map<unsigned int,unsigned int> mcl = tableC(observationMemberships);
      for(auto it = mcl.cbegin(); it != mcl.cend(); ++it)
        Rcpp::Rcout <<" - Size of cluster #"<< it->first << ": " << it->second << std::endl;
    }

    if (m_UseFence)
    {
      if (m_UseVerbose)
        Rcpp::Rcout << "Fence algorithm: "<< std::endl;

      iterativeFence(
        warpingParameters,
        iter,
        observationMemberships,
        observationDistances,
        m_WarpingPointer,
        m_OptimizerPointer,
        templateValues,
        warpedGrids,
        m_InputValues,
        templateGrids,
        m_DissimilarityPointer,
        clusterIndices,
        m_UseVerbose
      );

      if (m_UseVerbose)
        Rcpp::Rcout << "Done" << std::endl;
    }

    // Normalization
    if (m_UseVerbose)
      Rcpp::Rcout << "Parameter normalization: ";

    m_WarpingPointer->Normalize(warpingParameters, clusterIndices, observationMemberships);

    if (m_UseVerbose)
      Rcpp::Rcout << "Done" << std::endl;

    // Store parameter values in container
    warpingParametersContainer.slice(iter - 1) = warpingParameters;

    // Update individual warped grids
    if (m_UseVerbose)
      Rcpp::Rcout << "Update individual warped grids: ";

    warpedGrids = m_WarpingPointer->ApplyWarping(warpedGrids, warpingParameters);

    if (m_UseVerbose)
      Rcpp::Rcout << "Done" << std::endl;

    timer.step( "fence/norm/update "+ std::to_string(iter) );

    // Compute new templates
    if (m_UseVerbose)
      Rcpp::Rcout << "Compute new templates: " << std::endl;

    templateGridsContainer.slice(iter - 1) = templateGrids;
    templateValuesContainer(iter - 1) = templateValues;

    templateGrids.set_size(numberOfClusters, m_NumberOfPoints);
    templateValues.set_size(numberOfClusters, m_NumberOfDimensions, m_NumberOfPoints);

    this->UpdateTemplates(
      warpedGrids,
      clusterIndices,
      observationMemberships,
      templateGrids,
      templateValues
    );

    //check total smilarity
    if (m_CheckTotalSimilarity)
    {
      if (m_UseVerbose)
        Rcpp::Rcout << "Check total similarity: ";

      double totalDissimilarity = arma::sum(observationDistances);
      double oldTotalDissimilarity = arma::sum(oldObservationDistances);

      // if total distance increased
      if (oldTotalDissimilarity < totalDissimilarity)
      {
        totalDissimilarityCondition = false;
        templateGrids = templateGridsContainer.slice(iter - 1);
        templateValues = templateValuesContainer(iter - 1);
        observationDistances = oldObservationDistances;
        observationMemberships = oldObservationMemberships;

        if (m_UseVerbose)
          Rcpp::Rcout << "Total dissimilarity did not decrease. ";
      }

      if (m_UseVerbose)
        Rcpp::Rcout << "Done" << std::endl;
    }

    timer.step( "newtemplates "+ std::to_string(iter) );
  }

  warpingParametersContainer.resize(m_NumberOfObservations, numberOfParameters, iter);
  templateGridsContainer.slice(iter) = templateGrids;
  templateValuesContainer(iter) = templateValues;

  if (m_UseVerbose)
    Rcpp::Rcout << "End while iterations" << std::endl;

  //
  // output
  //
  if (m_UseVerbose)
    Rcpp::Rcout << "Final warping: ";

  arma::mat finalWarpingParameters = m_WarpingPointer->GetFinalWarping(
    warpingParametersContainer,
    observationMemberships,
    clusterIndices
  );

  if (m_UseVerbose)
    Rcpp::Rcout << "Done" << std::endl;

  // Convert cube to field for conversion to List in R
  arma::field<arma::mat> listOfEstimatedParameters(iter);
  arma::field<arma::mat> listOfTemplateGrids(iter);

  for (unsigned int k = 0;k < iter;++k)
  {
    listOfEstimatedParameters(k) = warpingParametersContainer.slice(k);
    listOfTemplateGrids(k) = templateGridsContainer.slice(k);
  }

  Rcpp::NumericVector out1 = Rcpp::wrap(overallCenter.centerGrid);
  out1.attr("dim") = R_NilValue;

  Rcpp::NumericVector out2 = Rcpp::wrap(overallCenter.distancesToCenter);
  out2.attr("dim") = R_NilValue;

  Rcpp::NumericVector out4 = Rcpp::wrap(observationDistances);
  out4.attr("dim") = R_NilValue;

  Rcpp::NumericVector out7 = Rcpp::wrap(observationMemberships + 1);
  out7.attr("dim") = R_NilValue;

  arma::field<arma::cube> out8 = templateValuesContainer.rows(0, iter);
  arma::field<arma::mat> out9 = listOfTemplateGrids;

  timer.step( "output ");

  return ListBuilder()
    .add("iterations",                  iter)
    .add("n_clust",                     m_NumberOfClusters)
    .add("overall_center_grid",         out1)
    .add("overall_center_values",       overallCenter.centerValues)
    .add("distances_to_overall_center", out2)
    .add("x_final",                     warpedGrids)
    .add("n_clust_final",               clusterIndices.size())
    .add("x_centers_final",             templateGrids)
    .add("y_centers_final",             templateValues)
    .add("templates_vec",               out8)
    .add("x_out_vec",                   out9)
    .add("labels",                      out7)
    .add("final_dissimilarity",         out4)
    .add("parameters_list",             listOfEstimatedParameters)
    .add("parameters",                  finalWarpingParameters)
    .add("timer",                       timer);
}

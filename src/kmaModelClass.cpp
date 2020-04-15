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
  Rcpp::Rcout << " - Initial seeds for cluster centers: " << m_SeedVector << std::endl;

  Rcpp::Rcout << "Information about the methods used within the algorithm:" << std::endl;
  Rcpp::Rcout << " - Warping method: " << warpingMethod << std::endl;
  Rcpp::Rcout << " - Center method: " << centerMethod << std::endl;
  Rcpp::Rcout << " - Dissimilarity method: " << dissimilarityMethod << std::endl;
  Rcpp::Rcout << " - Optimization method: " << optimizerMethod << std::endl;
  Rcpp::Rcout << " - Interpolation method: " << m_InterpolationMethod << std::endl;

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

void KmaModel::UpdateTemplates(const arma::mat& x_reg,
                               const arma::urowvec& ict,
                               const arma::urowvec& labels,
                               arma::cube& templates)
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

    for (unsigned int i = 0;i < ict.size();++i)
    {
      selectedObservations = arma::find(labels == ict(i));

      centerComputer = m_CenterPointer->GetCenter(
        x_reg.rows(selectedObservations),
        GetObservations(m_InputValues, selectedObservations),
        m_DissimilarityPointer
      );

      templates.tube(arma::span(i), arma::span::all) = centerComputer.centerValues;
    }
    break;

  case DistanceLoop:

    for (unsigned int i = 0;i < ict.size();++i)
    {
      selectedObservations = arma::find(labels == ict(i));

      centerComputer = m_CenterPointer->GetCenter(
        x_reg.rows(selectedObservations),
        GetObservations(m_InputValues, selectedObservations),
        m_DissimilarityPointer,
        m_NumberOfThreads
      );

      templates.tube(arma::span(i), arma::span::all) = centerComputer.centerValues;

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

  while (distanceCondition && membershipCondition && iterationCondition && totalDissimilarityCondition)
  {
    ++iter;

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

#ifdef _OPENMP
#pragma omp parallel for num_threads(m_NumberOfThreads)
#endif

    // inizializzo container warp_temp
    unsigned int numberOfTemplates = templates.n_rows;
    arma::rowvec index_temp(numberOfTemplates);
    arma::mat parameters_temp(numberOfParameters, numberOfTemplates);
    arma::colvec arg(numberOfParameters);
    arma::mat y_reg;
    arma::mat t_in;
    WarpingSet warpingSet;

    for (unsigned int i = 0;i < m_NumberOfObservations;++i)
    {
      y_reg = approx(
        x_reg.row(i),
        GetObservation(m_InputValues, i),
        m_InterpolationMethod
      );

      // Compute warping parameters for each template
      for (unsigned int j = 0;j < numberOfTemplates;++j)
      {
        t_in = templates(arma::span(j), arma::span::all, arma::span::all);

        if (m_NumberOfDimensions > 1)
          t_in = t_in.t();

        warpingSet = m_WarpingPointer->SetInputData(x_out, x_out, y_reg, t_in, m_DissimilarityPointer);

        auto fun = [this, &warpingSet] (const arma::vec &arg)
        {
          return this->m_WarpingPointer->GetDissimilarityAfterWarping(warpingSet, arg);
        };

        index_temp(j) = m_OptimizerPointer->Optimize(arg, m_WarpingPointer, fun);
        parameters_temp.col(j) = arg;
      }

      //fine iterazioni per ogni tempalte
      index(i) = index_temp.min();
      labels(i) = ict(arma::index_min(index_temp));
      parameters.col(i) = parameters_temp.col(arma::index_min(index_temp));
    }

    if (m_UseVerbose)
      Rcpp::Rcout << "Done" << std::endl;

    timer.step( "warping "+ std::to_string(iter) );

    //update current template list
    ict = unique(labels);

    //PRINT
    if (m_UseVerbose)
    {
      Rcpp::Rcout << "current cluster vector updated" << std::endl;
      ict.print();
      std::map<unsigned int,unsigned int> mcl = tableC(labels);
      for(auto it = mcl.cbegin(); it != mcl.cend(); ++it)
        Rcpp::Rcout <<"cluster num: "<< it->first << " has " << it->second << " elements;" << std::endl;
    }

    if (m_UseFence)
    {
      if (m_UseVerbose)
        Rcpp::Rcout << "Fence algorithm: "<< std::endl;

      iterativeFence(parameters, iter, labels, index, m_WarpingPointer, m_OptimizerPointer, templates,
                     x_reg, m_InputValues, x_out, m_DissimilarityPointer, ict, m_UseVerbose);

      if (m_UseVerbose)
        Rcpp::Rcout << "Done" << std::endl;
    }

    // normalizzazione
    if (m_UseVerbose)
      Rcpp::Rcout << "Parameter normalization: ";

    m_WarpingPointer->Normalize(parameters, ict, labels);

    if (m_UseVerbose)
      Rcpp::Rcout << "Done" << std::endl;

    // salvo parametri
    parameters_vec(arma::span::all, arma::span::all, arma::span(iter - 1)) = parameters;

    //update x_reg and x_out
    if (m_UseVerbose)
      Rcpp::Rcout << "Update x_reg and x_out: ";

    x_reg = m_WarpingPointer->ApplyWarping(x_reg, parameters);
    x_out = arma::linspace<arma::rowvec>(
      GetCommonLowerBound(x_reg),
      GetCommonUpperBound(x_reg),
      m_NumberOfPoints
    );

    if (m_UseVerbose)
      Rcpp::Rcout << "Done" << std::endl;

    x_out_vec(iter) = x_out;

    timer.step( "fence/norm/update "+ std::to_string(iter) );

    //compute new templates
    if (m_UseVerbose)
      Rcpp::Rcout << "Compute new templates: " << std::endl;

    templates_vec(iter - 1) = templates;
    templates.set_size(ict.size(), m_NumberOfPoints, m_NumberOfDimensions);

    this->UpdateTemplates(
      x_reg,
      ict,
      labels,
      templates
    );

    if (m_UseVerbose)
    {
      Rcpp::Rcout <<"Templates updated" << std::endl;
      Rcpp::Rcout << "While condition" << std::endl;
      Rcpp::Rcout << "dissim cambiata piu di toll: " << (sum( abs(index-index_old) < m_Tolerance) < m_NumberOfObservations) << std::endl;
      Rcpp::Rcout << "almeno un etichetta cambiata: " << (sum( labels == labels_old  ) != m_NumberOfObservations) << std::endl;
    }

    //check total smilarity
    if (m_CheckTotalSimilarity)
    {
      if (m_UseVerbose)
        Rcpp::Rcout << "Check total similarity: ";

      double tot = sum(index);
      double tot_old = sum(index_old);

      // sel la distanza totale aumenta
      if (tot_old < tot)
      {
        still_in = false;
        templates = templates_vec(iter - 1);
        index = index_old;
        labels = labels_old;
        x_out = x_out_vec(iter - 1);

        if (m_UseVerbose)
          Rcpp::Rcout << "Total similarity didn't increase. ";
      }

      if (m_UseVerbose)
        Rcpp::Rcout << "Done" << std::endl;
    }

    timer.step( "newtemplates "+ std::to_string(iter) );

  }//fine while

  parameters_vec.resize(numberOfParameters, m_NumberOfObservations, iter);
  templates_vec(iter) = templates;

  if (m_UseVerbose)
    Rcpp::Rcout << "End while iterations" << std::endl;

  //
  // output
  //
  if (m_UseVerbose)
    Rcpp::Rcout << "Final warping: ";

  arma::mat final_par = m_WarpingPointer->GetFinalWarping(parameters_vec, labels, ict);

  if (m_UseVerbose)
    Rcpp::Rcout << "Done" << std::endl;

  arma::field<arma::mat> par_vec(iter);
  for(unsigned int k = 0;k < iter;++k)
    par_vec(k) = parameters_vec.slice(k);

  Rcpp::NumericVector out1 = Rcpp::wrap(original_center.centerGrid);
  out1.attr("dim") = R_NilValue;

  Rcpp::NumericVector out2 = Rcpp::wrap(original_center.distancesToCenter);
  out2.attr("dim") = R_NilValue;

  Rcpp::NumericVector out3 = Rcpp::wrap(x_out);
  out3.attr("dim") = R_NilValue;

  Rcpp::NumericVector out4 = Rcpp::wrap(index);
  out4.attr("dim") = R_NilValue;

  Rcpp::NumericVector out7 = Rcpp::wrap(labels+1);
  out7.attr("dim") = R_NilValue;

  arma::field<arma::cube> out8 = templates_vec.cols(0, iter);
  arma::field<arma::rowvec> out9 = x_out_vec.cols(0, iter);

  timer.step( "output ");

  return ListBuilder()
    .add("iterations", iter)
    .add("n.clust", m_NumberOfClusters)
    .add("x.center.orig",out1)
    .add("y.center.orig",original_center.centerValues)
    .add("similarity.origin",out2)
    .add("x.final", x_reg)
    .add("n.clust.final", ict.size())
    .add("x.centers.final", out3)
    .add("y.centers.final",templates)
    .add("templates_vec",out8)
    .add("x_out_vec",out9)
    .add("labels",out7)
    .add("similarity.final",out4)
    .add("parameters.list", par_vec)
    .add("parameters", final_par)
    .add("timer",timer);
}

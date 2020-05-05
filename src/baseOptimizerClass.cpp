#include "baseOptimizerClass.h"

double BaseOptimizerFunction::AlignToTemplateCostFunction(unsigned n,
                                                          const double *x,
                                                          double *grad,
                                                          void *data)
{
  AlignToTemplateData *d = (AlignToTemplateData *) data;

  arma::rowvec params(n);
  for (unsigned int i = 0;i < n;++i)
    params(i) = x[i];

  WarpingSet warpingSet;
  warpingSet.inputGrid1 = d->inputGrid;
  warpingSet.inputValues1 = d->inputValues;
  warpingSet.inputGrid2 = d->templateGrid;
  warpingSet.inputValues2 = d->templateValues;
  warpingSet.dissimilarityPointer = d->dissimilarityPointer;

  double distanceValue = d->warpingPointer->GetDissimilarityAfterWarping(warpingSet, params);

  double resValue = (1.0 - d->penalizationWeight) * distanceValue * distanceValue;

  double normValue = arma::norm(params - d->warpingPointer->GetInitialPoint());
  resValue += d->penalizationWeight * normValue * normValue;

  return resValue;
}

double BaseOptimizerFunction::CenterTemplateCostFunction(unsigned n,
                                                         const double *x,
                                                         double *grad,
                                                         void *data)
{
  CenterTemplateData *d = (CenterTemplateData *) data;

  arma::rowvec params(n);
  for (unsigned int i = 0;i < n;++i)
    params(i) = x[i];

  WarpingSet warpingSet;
  warpingSet.dissimilarityPointer = d->dissimilarityPointer;
  warpingSet.inputGrid1 = d->templateGrid;
  warpingSet.inputValues1 = d->templateValues;

  unsigned int numberOfObservations = d->inputGrids.n_rows;
  double resValue = 0.0;

  for (unsigned int i = 0;i < numberOfObservations;++i)
  {
    warpingSet.inputGrid2 = d->inputGrids.row(i);
    warpingSet.inputValues2 = d->inputValues.row(i);
    double distanceValue = d->warpingPointer->GetDissimilarityAfterWarping(warpingSet, params);
    resValue += distanceValue * distanceValue;
  }

  resValue /= (double)numberOfObservations;
  resValue *= (1.0 - d->penalizationWeight);

  double normValue = arma::norm(params - d->warpingPointer->GetInitialPoint());
  resValue += d->penalizationWeight * normValue * normValue;

  return resValue;
}

double BaseOptimizerFunction::AlignToTemplate(arma::rowvec &initialParameters,
                                              const arma::rowvec &inputGrid,
                                              const arma::mat &inputValues,
                                              const arma::rowvec &templateGrid,
                                              const arma::mat &templateValues,
                                              const std::shared_ptr<BaseDissimilarityFunction> &dissimilarityPointer,
                                              const std::shared_ptr<BaseWarpingFunction> &warpingPointer)
{
    unsigned int numberOfParameters = warpingPointer->GetNumberOfParameters();
    nlopt_opt optimizer = this->GetOptimizer(numberOfParameters);

    arma::rowvec lowerBounds = warpingPointer->GetParameterLowerBounds();
    arma::rowvec upperBounds = warpingPointer->GetParameterUpperBounds();
    initialParameters = warpingPointer->GetInitialPoint();

    AlignToTemplateData extraData;
    extraData.inputGrid = inputGrid;
    extraData.inputValues = inputValues;
    extraData.templateGrid = templateGrid;
    extraData.templateValues = templateValues;
    extraData.dissimilarityPointer = dissimilarityPointer;
    extraData.warpingPointer = warpingPointer;
    extraData.penalizationWeight = m_PenalizationWeight;

    if (initialParameters.size() == 0)
        return this->AlignToTemplateCostFunction(
            numberOfParameters,
            &(initialParameters[0]),
            NULL,
            &extraData
        );

    nlopt_set_lower_bounds(optimizer, &(lowerBounds[0]));
    nlopt_set_upper_bounds(optimizer, &(upperBounds[0]));

    nlopt_set_min_objective(optimizer, this->AlignToTemplateCostFunction, &extraData);
    nlopt_set_xtol_rel(optimizer, m_ParameterRelativeTolerance);

    double fVal;
    int exitCode = nlopt_optimize(optimizer, &(initialParameters[0]), &fVal);

    nlopt_destroy(optimizer);

    if (exitCode < 0)
    {
      Rcpp::Rcout << fVal << " " << initialParameters << " " << lowerBounds << " " << upperBounds << std::endl;
      Rcpp::stop("NLOPT optimization failed.");
    }

    return fVal;
}

double BaseOptimizerFunction::CenterTemplate(arma::rowvec &initialParameters,
                                             const arma::rowvec &templateGrid,
                                             const arma::mat &templateValues,
                                             const arma::mat &inputGrids,
                                             const arma::cube &inputValues,
                                             const std::shared_ptr<BaseDissimilarityFunction> &dissimilarityPointer,
                                             const std::shared_ptr<BaseWarpingFunction> &warpingPointer)
{
  unsigned int numberOfParameters = warpingPointer->GetNumberOfParameters();
  nlopt_opt optimizer = this->GetOptimizer(numberOfParameters);

  arma::rowvec lowerBounds = warpingPointer->GetParameterLowerBounds();
  arma::rowvec upperBounds = warpingPointer->GetParameterUpperBounds();
  initialParameters = warpingPointer->GetInitialPoint();

  CenterTemplateData extraData;
  extraData.templateGrid = templateGrid;
  extraData.templateValues = templateValues;
  extraData.inputGrids = inputGrids;
  extraData.inputValues = inputValues;
  extraData.dissimilarityPointer = dissimilarityPointer;
  extraData.warpingPointer = warpingPointer;
  extraData.penalizationWeight = m_PenalizationWeight;

  if (initialParameters.size() == 0)
    return this->CenterTemplateCostFunction(
        numberOfParameters,
        &(initialParameters[0]),
        NULL,
        &extraData
    );

  nlopt_set_lower_bounds(optimizer, &(lowerBounds[0]));
  nlopt_set_upper_bounds(optimizer, &(upperBounds[0]));

  nlopt_set_min_objective(optimizer, this->CenterTemplateCostFunction, &extraData);
  nlopt_set_xtol_rel(optimizer, m_ParameterRelativeTolerance);

  double fVal;
  int exitCode = nlopt_optimize(optimizer, &(initialParameters[0]), &fVal);

  nlopt_destroy(optimizer);

  if (exitCode < 0)
    Rcpp::stop("NLOPT optimization failed.");

  return fVal;
}

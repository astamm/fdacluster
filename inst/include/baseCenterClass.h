#ifndef BASECENTERCLASS_H
#define BASECENTERCLASS_H

#include "baseDissimilarityClass.h"

#include <memory>

struct CenterType
{
    arma::rowvec centerGrid;
    arma::mat centerValues;
    arma::rowvec distancesToCenter;
};

class BaseCenterMethod
{
public:
  BaseCenterMethod()
  {
    m_SpanValue = 0.1;
    m_PolynomialDegree = 4;
  }
  virtual ~BaseCenterMethod() {};

  /// Compute center method.
  /**
   *  @param[inputGrid] Input grid on which observed functions are evaluated;
   *  @param[inputValues] Input function values on input grid;
   *  @param[distanceObject] Shared pointer to the base class Dissimilarity;
   *  @param[nbThreads] Number of threads to use during the computation.
   *
   *  @return A center object.
   */
  virtual CenterType GetCenter(
      const arma::mat& inputGrid,
      const arma::cube& inputValues,
      const std::shared_ptr<BaseDissimilarityFunction>& dissimilarityPointer
  ) = 0;

  virtual CenterType GetCenter(const arma::mat& inputGrid,
                               const arma::cube& inputValues,
                               const std::shared_ptr<BaseDissimilarityFunction>& dissimilarityPointer,
                               unsigned int nbThreads)
  {
    CenterType outputCenter;
    return outputCenter;
  }

  void SetSpanValue(const double &val) {m_SpanValue = val;}
  double GetSpanValue() {return m_SpanValue;}

  void SetPolynomialDegree(const unsigned int &val) {m_PolynomialDegree = val;}
  unsigned int GetPolynomialDegree() {return m_PolynomialDegree;}

private:
  double m_SpanValue;
  unsigned int m_PolynomialDegree;
};

#endif /* BASECENTERCLASS_H */

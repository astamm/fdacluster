#ifndef BASECENTERCLASS_H
#define BASECENTERCLASS_H

#include "dissimilarity.h"

#include <RcppArmadillo.h>

enum SpaceType
{
    Euclidean,
    UnitQuaternion
};

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
        m_SpanValue = 0.2;
        m_DeltaValue = 0.0;
        m_Space = Euclidean;
    }

    virtual ~BaseCenterMethod() {};

    /// Compute center method.
    /**
     *  @param[inputGrid] Input grid on which observed functions are evaluated;
     *  @param[inputValues] Input function values on input grid;
     *  @param[distanceObject] Shared pointer to the base class Dissimilarity;
     *
     *  @return A center object.
     */
    virtual CenterType GetCenter(
            const arma::mat& inputGrid,
            const arma::cube& inputValues,
            const std::shared_ptr<Dissimilarity>& dissimilarityPointer
    ) = 0;

    /// Compute center method in parallel (used if type of parallelization is 1).
    /**
     *  @param[inputGrid] Input grid on which observed functions are evaluated;
     *  @param[inputValues] Input function values on input grid;
     *  @param[distanceObject] Shared pointer to the base class Dissimilarity;
     *  @param[nbThreads] Number of threads to use during the computation.
     *
     *  @return A center object.
     */
    virtual CenterType GetCenterParallel(const arma::mat& inputGrid,
                                         const arma::cube& inputValues,
                                         const std::shared_ptr<Dissimilarity>& dissimilarityPointer,
                                         unsigned int nbThreads)
    {
        CenterType outputCenter;
        return outputCenter;
    }

    void SetSpanValue(const double num) {m_SpanValue = num;}
    double GetSpanValue() {return m_SpanValue;}

    void SetDeltaValue(const double num) {m_DeltaValue = num;};
    double GetDeltaValue() {return m_DeltaValue;}

    void SetSpace(const enum SpaceType &val) {m_Space = val;}
    enum SpaceType GetSpace() {return m_Space;}

private:
    double m_SpanValue;
    double m_DeltaValue;
    SpaceType m_Space;
};

#endif /* BASECENTERCLASS_H */

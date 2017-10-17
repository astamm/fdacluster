#ifndef CENTER_METHODS_HPP_
#define CENTER_METHODS_HPP_

#include<RcppArmadillo.h>
#include<memory>

#include "dissimilarity.h"
#include "utilities.h"

using namespace arma;

/// Center element returned by computeCenter and computeParallelCenter methods
struct center
{
    rowvec x_center;
    mat y_center;
    rowvec dissim_whit_origin;
};


/// Base class for all the center methods available.
class CenterMethod
{

public:
    /// Compute center method.
    /**
     *  @param[x] abscissa of the input functions;
     *  @param[y] values of the input functions;
     *  @param[dissim] shared pointer to the base class Dissimilarity.
     *  @param[x_out] grid of the returning center
     *
     *  @return A center element
     */
    virtual center computeCenter(const mat& x,const cube& y, std::shared_ptr<Dissimilarity>& dissim, const rowvec& x_out)=0;
    /// Compute center method in parallel (used if type of parallelization is 1).
    /**
     *  @param[x] abscissa of the input functions;
     *  @param[y] values of the input functions;
     *  @param[dissim] shared pointer to the base class Dissimilarity.
     *  @param[x_out] grid of the returning center
     *  @param[n_th] number of threads to use during the computation.
     *
     *  @return A center element
     */
    virtual center computeParallelCenter(const mat& x, const cube& y, std::shared_ptr<Dissimilarity>& dissim, const rowvec& x_out,uword n_th)=0;
    /// Set the parameters needed by the chosen center method.
    /**
     * @param[par] a row vector with options.
     */
    virtual void setParameters(const rowvec& par)=0;
};


/// Medoid center method finds the real medoid center.
class Medoid final: public CenterMethod
{
public:
    virtual center computeCenter(const mat& x, const cube& y, std::shared_ptr<Dissimilarity>& dissim, const rowvec& x_out);
    virtual center computeParallelCenter(const mat& x, const cube& y, std::shared_ptr<Dissimilarity>& dissim, const rowvec& x_out,uword n_th);
    virtual void setParameters(const rowvec& par) {};
};


/// PseudoMedoid center method finds the real medoid of each component and create a center.
class PseudoMedoid final: public CenterMethod
{
public:
    virtual center computeCenter(const mat& x, const cube& y, std::shared_ptr<Dissimilarity>& dissim, const rowvec& x_out);
    virtual center computeParallelCenter(const mat& x, const cube& y, std::shared_ptr<Dissimilarity>& dissim, const rowvec& x_out,uword n_th) {};
    virtual void setParameters(const rowvec& par) {};
};

/// Mean center method compute an approximation by lowess.
class Mean final: public CenterMethod
{
private:
    double span=0;
    double d=0;
public:
    virtual center computeCenter(const mat& x, const cube& y, std::shared_ptr<Dissimilarity>& dissim, const rowvec& x_out);
    virtual center computeParallelCenter(const mat& x, const cube& y, std::shared_ptr<Dissimilarity>& dissim, const rowvec& x_out,uword n_th){};
    virtual void setParameters(const rowvec& par);
};

#endif

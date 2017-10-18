#ifndef _WARPING_HPP
#define _WARPING_HPP

#include <RcppArmadillo.h>
#include <memory>
#include "dissimilarity.h"
#include "utilities.h"

/// Warping setting
/**
  * Input for warp member function.
  * It contains abscisa and values of the functions to warp and it is returned by
  * the set_function member function.
  */
struct warping_set
{
    rowvec xf;
    rowvec xg;
    mat yf;
    mat yg;
};


/// Base class for all The Warping Functions
/**
 * From this class all the warping available are derived.
 */
class WarpingFunction
{
protected:

    /// Pointer to dissimilarity object to use.
    std::shared_ptr<Dissimilarity> diss;

    /// Bounds for the warping.
    rowvec upperBound;
    rowvec lowerBound;

public:
    /// Member to create warping_set.
    /**
     * @param[x_f] abscissa f function;
     * @param[x_g] abscissa g function;
     * @param[y_f] values f function;
     * @param[y_g] values g function;
     * @param[d] pointer tho the base class Dissimilarity.
     *
     * @return a warping_set object.
     */
    warping_set set_function(const rowvec& x_f, const rowvec& x_g,
                             const mat& y_f, const mat& y_g,
                             std::shared_ptr<Dissimilarity>& d)
    {
        warping_set out;
        out.xf = x_f;
        out.xg = x_g;
        out.yf = y_f;
        out.yg = y_g;
        diss = d;
        return out;
    };

    /// Return upper bounds.
    rowvec get_upper_bound()
    {
        return upperBound;
    };
    /// Return lower bounds.
    rowvec get_lower_bound()
    {
        return lowerBound;
    };

    /// Apply warping to a matrix.
    /**
     * @param[x] abscissa to warp;
     * @param[par] warping parameters to apply;
     *
     * return abscissas warped;
     */
    virtual mat apply_warping(const mat& x, const mat& par )=0;

    /// Return number of parameters.
    virtual uword n_pars()=0;

    /// Set bounds given the input option different for each warping function.
    /**
     * @param[war_opt] input warping option.
     * @param[x] absissa to warp.
     */
    virtual void set_bounds(const rowvec& war_opt, const mat& x)=0;

    /// Set buonds given in a matrix.
    /**
     * @param[bou] bounds already computed;
     */
    void set_bounds(const mat& bou)
    {
        lowerBound = bou.row(0);
        upperBound = bou.row(1);
    }

    /// Compute final warping.
    /**
     * @param[parameters_vec] warping's parameters of each iteration;
     * @param[labels] final labels;
     * @param[ict] index current clusters;
     *
     * @return a matrix with the total warping parameters applied.
     */
    virtual mat final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict)=0;

    /// Normalize the warping parameters computed by clusters.
    /**
     * @param[par] warping parameters computed;
     * @param[ict] index current clusters;
     * @param[labels] current labels;
     */
    virtual void normalize(mat& par,const urowvec& ict,const urowvec& labels)=0;

    /// Compute dissimilarity after warp for optimization.
    /**
     * @param[w_set] warping_set element with the functions to warp.
     * @param[arg] best parameters that will be computed.
     */
    virtual double warp(const warping_set& w_set, const colvec& arg) const =0;
};


/// Affine warping function.
/**
 * An affine trasformation of the abscissa x given
 *  dilation and shift (d,s) is: d * x + s .
 */
class AffineFunction : public WarpingFunction
{

public:
    virtual uword n_pars()
    {
        return 2;
    }

    virtual mat apply_warping(const mat& x, const mat& par )
    {
        mat out(x.n_rows,x.n_cols);
        for(size_t obs=0; obs < x.n_rows; obs++)
            {
                out.row(obs) = par(0,obs)*x.row(obs) + par(1,obs);
            }
        return out;
    }
    virtual void set_bounds(const rowvec& war_opt, const mat& x)
    {
        double min_temp= min( util::uppers(x)- util::lowers(x) );
        double dl= war_opt(0), sl = war_opt(1);
        upperBound = { 1 + dl,  sl * min_temp};
        lowerBound = { 1 - dl, -sl * min_temp};
    }
    virtual mat final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict);

    virtual void normalize(mat& par,const urowvec& ict,const urowvec& labels);

    virtual double warp(const warping_set& w_set, const colvec& arg) const
    {
        return diss->compute(arg(0)*w_set.xf+arg(1),w_set.xg,w_set.yf,w_set.yg);
    }

};


/// Shift warping function.
/**
 * A shift trasformation of the abscissa x given
 * the shift ( s ) is: x + s .
 */
class ShiftFunction : public WarpingFunction
{
public:
    virtual uword n_pars()
    {
        return 1;
    }
    virtual mat apply_warping(const mat& x, const mat& par )
    {
        mat out(x.n_rows,x.n_cols);
        for(size_t obs=0; obs < x.n_rows; obs++)
            {
                out.row(obs) = x.row(obs) + par(0,obs);
            }
        return out;
    }

    virtual void set_bounds(const rowvec& war_opt, const mat& x)
    {
        double min_temp= min( util::uppers(x)- util::lowers(x) );
        double sl= war_opt[0];
        upperBound = { sl * min_temp };
        lowerBound = {-sl * min_temp };
    }
    virtual mat final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict);
    virtual void normalize(mat& par,const urowvec& ict,const urowvec& labels);

    virtual double warp(const warping_set& w_set, const colvec& arg) const
    {
        return diss->compute(w_set.xf+arg(0),w_set.xg,w_set.yf,w_set.yg);
    }

};


/// Dilation warping function.
/**
 * A dilation trasformation of the abscissa x given
 * the dialtion ( d ) is: d * x.
 */
class DilationFunction : public WarpingFunction
{
public:

    virtual uword n_pars()
    {
        return 1;
    }
    virtual mat apply_warping(const mat& x, const mat& par )
    {
        mat out(x.n_rows,x.n_cols);
        for(size_t obs=0; obs < x.n_rows; obs++)
            {
                out.row(obs) = par(0,obs)*x.row(obs);
            }
        return out;
    }
    virtual void set_bounds(const rowvec& war_opt, const mat& x)
    {
        double min_temp= min( util::uppers(x)- util::lowers(x) );
        double dl= war_opt[0];
        upperBound = { 1 + dl };
        lowerBound = { 1 - dl };
    }
    virtual mat final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict);
    virtual void normalize(mat& par,const urowvec& ict,const urowvec& labels);
    virtual double warp(const warping_set& w_set, const colvec& arg) const
    {
        return diss->compute(arg(0)*w_set.xf, w_set.xg,w_set.yf,w_set.yg);
    }
};


/// No warping function.
/**
 * The trasformation applied is the identity:
 */
class NoAlignmentFunction : public WarpingFunction
{
public:

    virtual uword n_pars()
    {
        return 0;
    }
    virtual mat apply_warping(const mat& x, const mat& par )
    {
        return x;
    }

    virtual void set_bounds(const rowvec& war_opt, const mat& x)
    {
        double min_temp= min( util::uppers(x)- util::lowers(x) );
        upperBound.set_size(0);
        lowerBound.set_size(0);
    }
    virtual mat final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict);
    virtual void normalize(mat& par,const urowvec& ict,const urowvec& labels) {};

    virtual double warp(const warping_set& w_set, const colvec& arg) const
    {
        return diss->compute(w_set.xf,w_set.xg,w_set.yf,w_set.yg);
    }
};

#endif

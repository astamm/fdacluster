// Copyright (C) 2017 Alessandro Zito (zito.ales@gmail.com)
//
// This file is part of Fdakmapp.
//
// Fdakmapp is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Fdakmapp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Fdakmapp.  If not, see <http://www.gnu.org/licenses/>.

#ifndef BASEWARPINGCLASS_H
#define BASEWARPINGCLASS_H

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
struct WarpingSet
{
  arma::rowvec xf;
  arma::rowvec xg;
  arma::mat yf;
  arma::mat yg;
};

/// Base class for all warping functions
/**
 * From this class all the warping available are derived.
 */
class BaseWarpingFunction
{
public:
  /// Member to create WarpingSet.
  /**
   * @param[x_f] abscissa f function;
   * @param[x_g] abscissa g function;
   * @param[y_f] values f function;
   * @param[y_g] values g function;
   *
   * @return a warping_set object.
   */
  WarpingSet SetInputData(
      const arma::rowvec &x_f,
      const arma::rowvec &x_g,
      const arma::mat &y_f,
      const arma::mat &y_g,
      const std::shared_ptr<Dissimilarity> &d
  );

  arma::rowvec GetParameterLowerBounds() {return m_ParameterLowerBounds;}
  arma::rowvec GetParameterUpperBounds() {return m_ParameterUpperBounds;}

  /// Apply warping to a matrix.
  /**
   * @param[x] abscissa to warp;
   * @param[par] warping parameters to apply;
   *
   * return abscissas warped;
   */
  virtual arma::mat apply_warping(const arma::mat &x, const arma::mat &par) = 0;

  /// Return number of parameters.
  virtual unsigned int GetNumberOfParameters() = 0;

  /// Set bounds given the input option different for each warping function.
  /**
   * @param[war_opt] input warping option.
   * @param[x] absissa to warp.
   */
  virtual void SetParameterBounds(const arma::rowvec &war_opt, const arma::mat &x) = 0;

  /// Set bounds given in a matrix.
  /**
   * @param[bou] bounds already computed;
   */
  void SetParameterBounds(const arma::mat &bou)
  {
    m_ParameterLowerBounds = bou.row(0);
    m_ParameterUpperBounds = bou.row(1);
  }

  /// Compute final warping.
  /**
   * @param[parameters_vec] warping's parameters of each iteration;
   * @param[labels] final labels;
   * @param[ict] index current clusters;
   *
   * @return a matrix with the total warping parameters applied.
   */
  virtual arma::mat GetFinalWarping(
      const arma::cube &parameters_vec,
      const arma::urowvec &labels,
      const arma::urowvec &ict) = 0;

  /// Normalize the warping parameters computed by clusters.
  /**
   * @param[par] warping parameters computed;
   * @param[ict] index current clusters;
   * @param[labels] current labels;
   */
  virtual void Normalize(
      arma::mat &par,
      const arma::urowvec &ict,
      const arma::urowvec &labels) = 0;

  /// Compute dissimilarity after warp for optimization.
  /**
   * @param[w_set] warping_set element with the functions to warp.
   * @param[arg] best parameters that will be computed.
   */
  virtual double GetDissimilarityAfterWarping(
      const WarpingSet &w_set,
      const arma::colvec &arg) = 0;

protected:
  /// Pointer to dissimilarity object to use.
  std::shared_ptr<Dissimilarity> m_DissimilarityPointer;

  /// Bounds for the warping.
  arma::rowvec m_ParameterLowerBounds;
  arma::rowvec m_ParameterUpperBounds;
};

/// Shift warping function.
/**
 * A shift transformation of the abscissa x given
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
    // double min_temp= min( util::uppers(x)- util::lowers(x) );
    double min_temp = arma::as_scalar(arma::min(arma::max(x, 1) - arma::min(x, 1)));
    double sl= war_opt[0];
    upperBound = { sl * min_temp };
    lowerBound = {-sl * min_temp };
  }
  virtual mat final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict);
  virtual void normalize(mat& par,const urowvec& ict,const urowvec& labels);

  virtual double warp(const warping_set& w_set, const colvec& arg) const
  {
    return diss->GetDistance(w_set.xf+arg(0),w_set.xg,w_set.yf,w_set.yg);
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
    double dl= war_opt[0];
    upperBound = { 1 + dl };
    lowerBound = { 1 - dl };
  }
  virtual mat final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict);
  virtual void normalize(mat& par,const urowvec& ict,const urowvec& labels);
  virtual double warp(const warping_set& w_set, const colvec& arg) const
  {
    return diss->GetDistance(arg(0)*w_set.xf, w_set.xg,w_set.yf,w_set.yg);
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
    upperBound.set_size(0);
    lowerBound.set_size(0);
  }
  virtual mat final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict);
  virtual void normalize(mat& par,const urowvec& ict,const urowvec& labels) {};

  virtual double warp(const warping_set& w_set, const colvec& arg) const
  {
    return diss->GetDistance(w_set.xf,w_set.xg,w_set.yf,w_set.yg);
  }
};

#endif /* BASEWARPINGCLASS_H */

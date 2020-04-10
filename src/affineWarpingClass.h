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

#ifndef AFFINEWARPINGCLASS_H
#define AFFINEWARPINGCLASS_H

#include <RcppArmadillo.h>
#include <memory>

#include "baseWarpingClass.h"
#include "dissimilarity.h"
#include "utilities.h"

class AffineWarpingFunction : public BaseWarpingFunction
{
  /**
   * An affine transformation of the abscissa x given
   *  dilation and shift (d,s) is: d * x + s .
   */

public:
  unsigned int GetNumberOfParameters() {return 2;}

  arma::mat ApplyWarping(const arma::mat &x, const arma::mat &par);
  void SetParameterBounds(const arma::rowvec &war_opt, const arma::mat &x);

  arma::mat GetFinalWarping(
      const arma::cube &parameters_vec,
      const arma::urowvec &labels,
      const arma::urowvec &ict
  );

  void Normalize(
      arma::mat &par,
      const arma::urowvec &ict,
      const arma::urowvec &labels
  );

  double GetDissimilarityAfterWarping(
      const WarpingSet &w_set,
      const arma::colvec &arg
  );
};

#endif /* AFFINEWARPINGCLASS_H */

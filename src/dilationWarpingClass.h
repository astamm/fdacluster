#ifndef DILATIONWARPINGCLASS_H
#define DILATIONWARPINGCLASS_H

#include "baseWarpingClass.h"

class DilationWarpingFunction : public BaseWarpingFunction
{
  /**
   * A dilation trasformation of the abscissa x given
   * the dialtion ( d ) is: d * x.
   */

public:
  unsigned int GetNumberOfParameters() {return 1;}

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
      const WarpingSet &warpingSet,
      const arma::colvec &arg
  );
};

#endif /* DILATIONWARPINGCLASS_H */

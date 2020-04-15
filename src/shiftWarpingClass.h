#ifndef SHIFTWARPINGCLASS_H
#define SHIFTWARPINGCLASS_H

#include "baseWarpingClass.h"

class ShiftWarpingFunction : public BaseWarpingFunction
{
  /**
   * A shift transformation of the abscissa x given
   * the shift ( s ) is: x + s .
   */

public:
  unsigned int GetNumberOfParameters();

  arma::mat ApplyWarping(const arma::mat &x, const arma::mat &par);
  void SetParameterBounds(const arma::rowvec &warpingOptions, const arma::mat &x);

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

#endif /* SHIFTWARPINGCLASS_H */

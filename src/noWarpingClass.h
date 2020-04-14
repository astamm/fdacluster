#ifndef NOWARPINGCLASS_H
#define NOWARPINGCLASS_H

#include "baseWarpingClass.h"

class NoWarpingFunction : public BaseWarpingFunction
{
  /**
   * The applied transformation is the identity:
   */

public:
  unsigned int GetNumberOfParameters() {return 0;}

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

#endif /* NOWARPINGCLASS_H */

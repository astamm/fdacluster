#ifndef DISSIMILARITY_HPP_
#define DISSIMILARITY_HPP_

#include<RcppArmadillo.h>

using namespace arma;

/// Element with approximated abiscissa and approximated functions to compute dissimilarity.
struct grid
{
    mat yf_sim;
    mat yg_sim;
    rowvec x_sim;
};


/// Base class for all the available dissimilarity.
class Dissimilarity
{
public:

    Dissimilarity() {};
    /// compute dissimilarity method different for each derived class
    virtual double compute(const rowvec& xf, const rowvec& xg,const mat& yf, const mat& yg)=0;
    /// compute common grid and approximations to compute dissimilarity
    grid setGrid(const rowvec& xf, const rowvec& xg, const mat& yf, const mat& yg);
};


/// Pearson similarity
/**
  *The function is returned negative to be a dissimilarity and
  * fit the minimization of the algorithm.
  */
 class Pearson final: public Dissimilarity
{
public:

    Pearson():Dissimilarity() {};

    double compute(const rowvec& xf, const rowvec& xg,
                   const mat& yf, const mat& yg);
};

/// L2 Distance
class L2 final: public Dissimilarity
{
public:
    L2():Dissimilarity() {};
    double compute(const rowvec& xf, const rowvec& xg,
                   const mat& yf, const mat& yg);
};



#endif

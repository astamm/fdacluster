#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_

#include <RcppArmadillo.h>

#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <iterator>


namespace util{

/// which_out Find the element out of the range.
/**
 *  \param vec vector of double.
 *  \param lo double lower bound.
 *  \param up  double upper bound.
 *
 *  \return vector of the position of the element out of the range.
 */

arma::urowvec which_out(arma::rowvec vec, double lo,double up);




/// quantile: Find the requested quantile.
/**
 *  \param vec vector of doubles.
 *  \param per  requested quantile (between 0 and 1).
 *  \return value in the requested position.
 */
double quantile(const arma::rowvec vec, const double per);



/// normalize NA in multidimensional function.
/**
 *  \param  mat  matrix representing a multidimensional function.
 *  \return matrix with NA normalized.
 */
arma::mat norm_ex(const arma::mat& y);


/// Created a ziiped vector.
/** Fill the zipped vector with pairs consisting of the
*   corresponding elements of a and b. (This assumes
*   that the vectors have equal length).
*/
template <typename A, typename B>
void zip(
    const std::vector<A> &a,
    const std::vector<B> &b,
    std::vector<std::pair<A,B>> &zipped)
{
    for(size_t i=0; i<a.size(); ++i)
    {
        zipped.push_back(std::make_pair(a[i], b[i]));
    }
}

/// Unzip a zipped vector.
/** Write the first and second element of the pairs in
* the given zipped vector into a and b. (This assumes
* that the vectors have equal length)
*/
template <typename A, typename B>
void unzip(
    const std::vector<std::pair<A, B>> &zipped,
    std::vector<A> &a,
    std::vector<B> &b)
{
    for(size_t i=0; i<a.size(); i++)
    {
        a[i] = zipped[i].first;
        b[i] = zipped[i].second;
    }
}



/// compute max values by row avoiding NA value
/**
 *  \param x matrix of doubles.
 *  \return vector of upper values.
 */
arma::rowvec uppers(const arma::mat& x);


/// compute min values by row avoiding NA value
/**
 *  \param x matrix of doubles.
 *  \return vector of lower values.
 */
arma::rowvec lowers(const arma::mat& x);


/// extract an obsevation from a cube
/**
 *  \param  y cube of observations.
 *  \parm i index of the observation to extract.
 *  \return etracted observation.
 */
const arma::mat observation(const arma::cube& y, arma::uword i);


/// extract many obsevations from a cube
/**
 *  \param y cube.
 *  \parm ind indeces of the observations to extract.
 *  \return extracted observations.
 */
const arma::cube observations(const arma::cube& y, arma::urowvec ind);


/// extract requested ebscissa.
/**
 *  \param mat  matrix of grids.
 *  \parm i index of the abscissa to extract.
 *  \return extracted abscissa.
 */
const arma::rowvec abscissa(const arma::mat& x,arma::uword i);

/// Approximate a function on a new grid.
/**
 *  \param x  original absicssa.
 *  \param y fuction to approximate.
 *  \return xx new abscissa.
 */
arma::mat approx(const arma::rowvec& x,
           const arma::mat& y,
           const arma::rowvec& xx);

/// R table function in c++
std::map<arma::uword, arma::uword> tableC(arma::urowvec x);


}

#endif

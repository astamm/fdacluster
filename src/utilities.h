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

#ifndef UTILITIES_HPP_
#define UTILITIES_HPP_

#include <RcppArmadillo.h>


#include<memory>
#include <vector>
#include <string>
#include <algorithm>
#include <iostream>
#include <iterator>


namespace util
{

/// which_out Find the element out of the range.
/**
 *  @param[vec] vector of double.
 *  @param[lo] double lower bound.
 *  @param[up]  double upper bound.
 *
 *  @return  A vector of the position of the element out of the range.
 */

arma::urowvec which_out(arma::rowvec vec, double lo,double up);




/// quantile: Find the requested quantile.
/**
 *  @param[vec] vector of doubles.
 *  @param[per]  requested quantile (between 0 and 1).
 *
 *  @return value in the requested position.
 */
double quantile(const arma::rowvec vec, const double per);



/// normalize NA in multidimensional function.
/**
 *  @param[y]  matrix representing a multidimensional function.
 *  @return matrix with NA normalized.
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
 *  @param[x] matrix of doubles.
 *  @return vector of upper values.
 */
arma::rowvec uppers(const arma::mat& x);


/// compute min values by row avoiding NA value
/**
 *  @param x matrix of doubles;
 *
 *  @return vector of lower values.
 */
arma::rowvec lowers(const arma::mat& x);


/// extract an obsevation from a cube
/**
 *  @param[y] cube of observations;
 *  @param[i] index of the observation to extract;
 *
 *  @return extracted observation.
 */
const arma::mat observation(const arma::cube& y, arma::uword i);


/// extract many obsevations from a cube
/**
 *  @param[y] cube.
 *  @parm[ind] indeces of the observations to extract.
 *  @return extracted observations.
 */
const arma::cube observations(const arma::cube& y, arma::urowvec ind);


/// extract requested ebscissa.
/**
 *  @param[mat]  matrix of grids.
 *  @param[i] index of the abscissa to extract.
 *  @return extracted abscissa.
 */
const arma::rowvec abscissa(const arma::mat& x,arma::uword i);

/// Approximate a function on a new grid.
/**
 *  @param[x]  original absicssa.
 *  @param[x] fuction to approximate.
 *
 *  @return xx new abscissa.
 */
arma::mat approx(const arma::rowvec& x,
                 const arma::mat& y,
                 const arma::rowvec& xx);

/// R table function in c++
/**
 * @param[x] input vector;
 *
 * @return the table of the input vector.
 */
std::map<arma::uword, arma::uword> tableC(arma::urowvec x);



/// List builder for build big list to return to R
class ListBuilder
{

public:

    ListBuilder() {};
    ~ListBuilder() {};

    inline ListBuilder& add(const std::string& name, SEXP x)
    {
        names.push_back(name);
        elements.push_back(PROTECT(x));
        return *this;
    }

    template <typename T>
    inline ListBuilder& add(const std::string& name, const T& x)
    {
        names.push_back(name);
        elements.push_back(PROTECT(Rcpp::wrap(x)));
        return *this;
    }

    inline operator Rcpp::List() const
    {
        Rcpp::List result(elements.size());
        for (size_t i = 0; i < elements.size(); ++i)
            {
                result[i] = elements[i];
            }
        result.attr("names") = Rcpp::wrap(names);
        UNPROTECT(elements.size());
        return result;
    }

    inline operator Rcpp::DataFrame() const
    {
        Rcpp::List result = static_cast<Rcpp::List>(*this);
        result.attr("class") = "data.frame";
        result.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, XLENGTH(elements[0]));
        return result;
    }

private:

    std::vector<std::string> names;
    std::vector<SEXP> elements;

    ListBuilder(ListBuilder const&) {};

};


/// Factory class
template<typename D>
class SharedFactory
{

public:
    typedef std::unordered_map< std::string, std::function< std::shared_ptr<D>() > > registry_map;

    registry_map map;

    // use this to instantiate the proper Derived class
    std::shared_ptr<D> instantiate(const std::string& name)
    {
        auto it = map.find(name);
        return it == map.end() ? nullptr : (it->second)();
    }

    template<typename T>
    void FactoryRegister(std::string name)
    {
        map[name] = []()
        {
            return std::make_shared<T>();
        };
        //std::cout << "Registering class '" << name << "'\n";
    }

};


}

#endif

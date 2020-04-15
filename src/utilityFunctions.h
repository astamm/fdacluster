#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include <RcppArmadillo.h>
#include <functional> // for std::function
#include <memory> // for std::shared_ptr...
#include <map>
#include <string>
#include <unordered_map> // for std::unordered_map
#include <utility>
#include <vector>

/// R table function in c++
/**
 * @param[inputLabels] Vector of input lables to count.
 *
 * @return Table containing the label counts.
 */
std::map<unsigned int, unsigned int> tableC(const arma::urowvec &inputLabels);

/// Find out-of-range elements.
/**
 *  @param[inputValues] Input values to be filtered.
 *  @param[lowerBound] Lower bound of included values.
 *  @param[upperBound] Upper bound of included values.
 *
 *  @return Indices of out-of-range elements.
 */
arma::urowvec which_out(const arma::rowvec &inputValues, const double lowerBound, const double upperBound);

/// Find the requested quantile.
/**
 *  @param[inputValues] Vector of input values.
 *  @param[quantileOrder] Order of the requested quantile (between 0 and 1).
 *
 *  @return Requested quantile.
 */
double quantile(const arma::rowvec &inputValues, const double quantileOrder);

/// Create a zipped vector.
/** Fill the zipped vector with pairs consisting of the
 *   corresponding elements of a and b. (This assumes
 *   that the vectors have equal length).
 */
template <typename A, typename B>
void zip(const std::vector<A> &a,
         const std::vector<B> &b,
         std::vector< std::pair<A,B> > &zipped)
{
    unsigned int vecDim = a.size();

    if (b.size() != a.size())
        Rcpp::stop("a and b should be of same size.");

    zipped.resize(vecDim);

    for (unsigned int i = 0;i < vecDim;++i)
        zipped[i] = std::make_pair(a[i], b[i]);
}

/// Unzip a zipped vector.
/** Write the first and second element of the pairs in
 * the given zipped vector into a and b. (This assumes
 * that the vectors have equal length)
 */
template <typename A, typename B>
void unzip(const std::vector< std::pair<A, B> > &zipped,
           std::vector<A> &a,
           std::vector<B> &b)
{
    unsigned int vecDim = zipped.size();
    a.resize(vecDim);
    b.resize(vecDim);

    for (unsigned i = 0;i < vecDim;++i)
    {
        a[i] = zipped[i].first;
        b[i] = zipped[i].second;
    }
}

/// Find the maximal lower bound of a set of grids.
/**
 *  @param[gridValues] Matrix of input grids.
 *
 *  @return Common lower bound.
 */
double GetCommonLowerBound(const arma::mat& gridValues);

/// Find the minimal upper bound of a set of grids.
/**
 *  @param[gridValues] Matrix of input grids.
 *
 *  @return Common upper bound.
 */
double GetCommonUpperBound(const arma::mat& gridValues);

/// Extract one observation from a cube.
/**
 *  @param[inputData] Data array in arma::cube format.
 *  @param[observationIndex] Index of the observation to extract.
 *
 *  @return Extracted observation.
 */
arma::mat GetObservation(const arma::cube& inputData, unsigned int observationIndex);

/// Extract several obsevations from a cube.
/**
 *  @param[inputData] Data array in arma::cube format.
 *  @param[observationIndices] Indices of the observations to extract.
 *  @return Extracted observations.
 */
arma::cube GetObservations(const arma::cube& inputData, arma::uvec& observationIndices);

/// Approximate a function on a new grid by linear interpolation.
/**
 *  @param[inputGrid] Original grid.
 *  @param[inputValues] Original values.
 *  @param[outputGrid] Interpolating grid.
 *
 *  @return Interpolated values.
 */
Rcpp::List approx(const arma::rowvec& inputGrid,
                  const arma::mat& inputValues,
                  const std::string interpolationMethod = "Linear"
);

/// List builder for build big list to return to R
class ListBuilder
{
public:
    ListBuilder() {};
    ~ListBuilder() {};

    inline ListBuilder& add(const std::string& name, SEXP x)
    {
        m_Names.push_back(name);
        m_Elements.push_back(PROTECT(x));

        return *this;
    }

    template <typename T>
    inline ListBuilder& add(const std::string& name, const T& x)
    {
        m_Names.push_back(name);
        m_Elements.push_back(PROTECT(Rcpp::wrap(x)));

        return *this;
    }

    inline operator Rcpp::List() const
    {
        Rcpp::List result(m_Elements.size());

        for (unsigned int i = 0;i < m_Elements.size();++i)
            result[i] = m_Elements[i];

        result.attr("names") = Rcpp::wrap(m_Names);
        UNPROTECT(m_Elements.size());

        return result;
    }

    inline operator Rcpp::DataFrame() const
    {
        Rcpp::List result = static_cast<Rcpp::List>(*this);
        result.attr("class") = "data.frame";
        result.attr("row.names") = Rcpp::IntegerVector::create(NA_INTEGER, XLENGTH(m_Elements[0]));
        return result;
    }

private:
    ListBuilder(ListBuilder const&) {};

    std::vector<std::string> m_Names;
    std::vector<SEXP> m_Elements;
};

/// Factory class
template <typename ObjectType>
class SharedFactory
{
public:
    using RegistryMap = std::unordered_map<std::string, std::function< std::shared_ptr<ObjectType>() > >;

    // use this to instantiate the proper Derived class
    std::shared_ptr<ObjectType> Instantiate(const std::string &name)
    {
        auto it = m_Map.find(name);
        return it == m_Map.end() ? nullptr : (it->second)();
    }

    template <typename T>
    void Register(std::string name)
    {
        m_Map[name] = []()
        {
            return std::make_shared<T>();
        };
    }

private:
    RegistryMap m_Map;
};

#endif /* UTILITYFUNCTIONS_H */

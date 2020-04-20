#ifndef UTILITYFUNCTIONS_H
#define UTILITYFUNCTIONS_H

#include <RcppArmadillo.h>
#include <functional> // for std::function
#include <memory> // for std::shared_ptr...
#include <map>
#include <string>
#include <unordered_map>
#include <utility>
#include <vector>

/// R table function in C++
/**
 * @param[inputLabels] Vector of input labels to count.
 *
 * @return Table containing the label counts.
 */
std::map<unsigned int, unsigned int> tableCpp(const arma::urowvec &inputLabels);

/// Extract several observations from a cube.
/**
 *  @param[inputData] Data array in arma::cube format.
 *  @param[observationIndices] Indices of the observations to extract.
 *  @return Extracted observations.
 */
arma::cube GetObservations(const arma::cube& inputData, arma::uvec& observationIndices);

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
    using RegistryMap = std::unordered_map<std::string, std::function<std::shared_ptr<ObjectType>()> >;

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

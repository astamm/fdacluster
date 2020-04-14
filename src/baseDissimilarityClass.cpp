#include "baseDissimilarityClass.h"

#include <squad.h>

FunctionPairType BaseDissimilarityFunction::GetComparableFunctions(const arma::rowvec& grid1,
                                                                   const arma::rowvec& grid2,
                                                                   const arma::mat& values1,
                                                                   const arma::mat& values2)
{
    FunctionPairType outputPair;
    outputPair.Grid.reset();
    outputPair.Values1.reset();
    outputPair.Values2.reset();

    // dimensions of inputs function
    unsigned int nDim = values1.n_rows;

    if (values2.n_rows != nDim)
        Rcpp::stop("Function domains have not the same dimension.");

    unsigned int nPts1 = values1.n_cols;

    if (grid1.size() != nPts1)
        Rcpp::stop("Function 1 is not properly evaluated on its grid.");

    unsigned int nPts2 = values2.n_cols;

    if (grid2.size() != nPts2)
        Rcpp::stop("Function 2 is not properly evaluated on its grid.");

    //
    // Missing value detection and elimination
    //

    arma::mat cleanValues1(nDim, nPts1, arma::fill::zeros), cleanValues2(nDim, nPts2, arma::fill::zeros);
    arma::rowvec cleanGrid1(nPts1, arma::fill::zeros), cleanGrid2(nPts2, arma::fill::zeros);

    unsigned int c = 0;
    for (unsigned int i = 0;i < nPts1;++i)
    {
        if (arma::is_finite(grid1(i)) && arma::is_finite(values1.col(i)))
        {
            cleanGrid1(c) = grid1(i);
            cleanValues1.col(c) = values1.col(i);
            ++c;
        }
    }

    cleanGrid1.resize(c);
    cleanValues1.resize(nDim, c);

    c = 0;
    for (unsigned int i = 0;i < nPts2;++i)
    {
        if (arma::is_finite(grid2(i)) && arma::is_finite(values2.col(i)))
        {
            cleanGrid2(c) = grid2(i);
            cleanValues2.col(c) = values2.col(i);
            ++c;
        }
    }

    cleanGrid2.resize(c);
    cleanValues2.resize(nDim, c);

    if (cleanValues1.n_cols == 0 || cleanValues2.n_cols == 0)
    {
        Rcpp::warning("Dissimilarity: at least one function only contains missing values.");
        return outputPair;
    }

    // Compute smallest common grid

    double xMin = std::max(cleanGrid1.min(), cleanGrid2.min());
    double xMax = std::min(cleanGrid1.max(), cleanGrid2.max());

    if (xMin >= xMax)
    {
        Rcpp::warning("Disimilarity: domain intersection is empty.");
        return outputPair;
    }

    unsigned int nPts = std::min(nPts1, nPts2);
    arma::rowvec xCommon = arma::linspace<arma::rowvec>(xMin, xMax, nPts);
    outputPair.Grid = xCommon;

    if (m_Space == Euclidean)
    {
        arma::rowvec workVector;

        for (unsigned int i = 0;i < nDim;++i)
        {
            arma::interp1(cleanGrid1, cleanValues1.row(i), xCommon, workVector, "*linear");
            outputPair.Values1.row(i) = workVector;
            arma::interp1(cleanGrid2, cleanValues2.row(i), xCommon, workVector, "*linear");
            outputPair.Values2.row(i) = workVector;
        }
    }
    else if (m_Space == UnitQuaternion)
    {
        outputPair.Values1 = Rcpp::as<arma::mat>(squad::RegularizeGrid(
            Rcpp::wrap(cleanGrid1),
            Rcpp::wrap(cleanValues1),
            xMin, xMax, nPts
        ));

        outputPair.Values2 = Rcpp::as<arma::mat>(squad::RegularizeGrid(
            Rcpp::wrap(cleanGrid2),
            Rcpp::wrap(cleanValues2),
            xMin, xMax, nPts
        ));
    }
    else
        Rcpp::Rcout << "Interpolation operations for the requested space are not yet implemented." << std::endl;


    return outputPair;
}

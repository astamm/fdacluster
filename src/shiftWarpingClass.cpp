#include "shiftWarpingClass.h"

unsigned int ShiftWarpingFunction::GetNumberOfParameters()
{
    return 1;
}

arma::mat ShiftWarpingFunction::ApplyWarping(const arma::mat &x, const arma::mat &par)
{
    arma::mat out(x.n_rows,x.n_cols);

    for (unsigned int i = 0;i < x.n_rows;++i)
        out.row(i) = x.row(i) + par(0, i);

    return out;
}

void ShiftWarpingFunction::SetParameterBounds(const arma::rowvec &war_opt,
                                              const arma::mat &x)
{
    double workScalar = arma::as_scalar(arma::min(arma::max(x, 1) - arma::min(x, 1)));
    double sl = war_opt(0);
    m_ParameterLowerBounds = { -sl * workScalar };
    m_ParameterUpperBounds = {  sl * workScalar };
}

arma::mat ShiftWarpingFunction::GetFinalWarping(const arma::cube &parameters_vec,
                                                const arma::urowvec &labels,
                                                const arma::urowvec &ict)
{
    unsigned int numberOfParameters = parameters_vec.n_rows;
    unsigned int numberOfObservations = parameters_vec.n_cols;
    unsigned int numberOfIterations = parameters_vec.n_slices;
    arma::mat out(numberOfParameters, numberOfObservations);
    out.row(0).zeros();
    arma::rowvec a;

    for (unsigned int i = 0;i < numberOfIterations;++i)
    {
        a = parameters_vec(arma::span(0), arma::span::all, arma::span(i));
        out.row(0) = out.row(0) + a;
    }

    arma::urowvec sel;
    arma::colvec m;

    for (unsigned int k = 0;k < ict.size();++k)
    {
        sel = arma::find(labels == ict(k)).t();
        m = arma::mean(out.cols(sel), 1);

        for (unsigned int i = 0;i < sel.size();++i)
            out(0, sel(i)) = out(0, sel(i)) - m(0);
    }

    return out;
}

void ShiftWarpingFunction::Normalize(arma::mat &par,
                                     const arma::urowvec &ict,
                                     const arma::urowvec &labels)
{
    arma::urowvec sel;
    arma::colvec par_mean;

    for (unsigned int i = 0;i < ict.size();++i)
    {
        sel = arma::find(labels == ict(i)).t();

        // calcolo medie cluster
        par_mean = arma::mean(par.cols(sel), 1);

        // aggiorno shift e dilation
        for (unsigned int j = 0;j < sel.size();++j)
            par(0, sel(j)) = par(0, sel(j)) - par_mean(0);
    }
}

double ShiftWarpingFunction::GetDissimilarityAfterWarping(const WarpingSet &warpingSet,
                                                          const arma::colvec &arg)
{
    return warpingSet.dissimilarityPointer->GetDistance(
            warpingSet.inputGrid1 + arg(0),
            warpingSet.inputGrid2,
            warpingSet.inputValues1,
            warpingSet.inputValues2
    );
}

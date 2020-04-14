#include "newCentersMethod.h"
#include "utilities.h"

#ifdef _OPENMP
#include <omp.h>
#endif

void newCenters(const arma::mat& x_reg,
                const arma::cube& y,
                const arma::rowvec& x_out,
                std::shared_ptr<BaseDissimilarityFunction>& dissim,
                std::shared_ptr<BaseCenterMethod>& cen,
                const arma::urowvec& par_opt,
                arma::cube& templates,
                const arma::urowvec& ict,
                const arma::urowvec& labels,
                const bool show_iter
)
{
    // switch to choose how to parallelize
    // case 0 trivial  each thread one cluster
    //case 1 each cluster all the threads (available only with medoid)

    switch(par_opt(1))
    {
    case 0:

#ifdef _OPENMP
#pragma omp parallel for num_threads(par_opt(0))
#endif

        for (unsigned int i = 0;i < ict.size();++i)
        {
            arma::urowvec sel = arma::find(labels == ict(i)).t();

            CenterType a = cen->GetCenter(
                x_reg.rows(sel),
                util::GetObservations(y, sel),
                dissim
            );

            templates.tube(arma::span(i), arma::span::all) = a.centerValues.t();
        }
        break;

    case 1:

        for (unsigned int i = 0;i < ict.size();++i)
        {
            arma::urowvec sel = arma::find(labels == ict(i)).t();

            CenterType a = cen->GetCenterParallel(
                x_reg.rows(sel),
                util::GetObservations(y, sel),
                dissim,
                par_opt(0)
            );

            templates.tube(arma::span(i), arma::span::all) = a.centerValues.t();

            if (show_iter)
                Rcpp::Rcout << "Template num. " << i << " computed" << std::endl;
        }
        break;
    }

}// fine new centers

#include "RcppArmadillo.h"

#include "factory.hpp"
#include "dissimilarity.hpp"
#include "warping.hpp"
#include "utilities.hpp"
#include "kma_model.hpp"
#include "checkin.hpp"

Rcpp::List kmap(const arma::mat& x,
                const arma::cube& y,
                const arma::rowvec& seeds,
                const arma::uword n_clust,
                const std::string warping_method,
                const std::string center_method,
                const std::string similarity_method,
                const std::string optim_method,
                const arma::rowvec& warping_opt,
                const arma::rowvec& center_opt,
                const arma::rowvec& out_opt,
                const bool fence,
                const bool check_total_similarity,
                const bool show_iter,
                const bool comp_original_center,
                const arma::urowvec par_opt)
{

    KmaModel modello(x, y, n_clust, warping_method, center_method, similarity_method, optim_method, seeds,
                     warping_opt, center_opt, out_opt(0),out_opt(1), fence, out_opt(2), show_iter, check_total_similarity, comp_original_center, par_opt);

    return  modello.execute();

}

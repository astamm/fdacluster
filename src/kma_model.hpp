#ifndef _KMA_MODEL_HPP
#define _KMA_MODEL_HPP

#include <RcppArmadillo.h>

#include "dissimilarity.hpp"
#include "warping.hpp"
#include "center_methods.hpp"
#include "optimizer.hpp"


using  namespace arma;

/// Main class.
/** This class handles loading of the problem and execution of the algorithm.
 */
class KmaModel
{

    const mat x;
    const cube y;
    const uword n_clust;
    const rowvec seeds;

    std::shared_ptr<Dissimilarity> dissim;
    std::shared_ptr<CenterMethod> cen;
    std::shared_ptr<WarpingFunction> warping;
    std::shared_ptr<OptimizerMethod> optimizer;

    const std::string optim_method;
    const rowvec warping_opt;
    const uword n_out;
    const double toll;
    const uword iter_max;
    bool fence;
    bool check_total_similarity;
    bool show_iter;
    bool com_oc;
    uword n_obs;
    uword n_camp;
    uword n_dim;
    const urowvec parallel_opt;

public:
    /// Constructor that load the problem.
    KmaModel(const mat& t_x, const cube& t_y,
             uword t_n_clust, std::string t_warping_method, std::string t_center_method,
             std::string t_similarity_method, std::string t_optim_method,
             const rowvec& t_seeds,
             const rowvec& t_warping_opt, const rowvec& t_center_opt, double t_n_out,
             double t_toll, bool t_fence, uword t_iter_max,bool t_show_iter,
             bool t_check_total_similarity,bool comp_original_center, urowvec par_opt);

    /// Method the execute the algorithm.
    Rcpp::List execute();

};

#endif

#include <RcppArmadillo.h>

#include "kma_model.h"
#include "checkin.h"

//
// KmaModle constructor
//

KmaModel::KmaModel(
    const mat& t_x, const cube& t_y,
    uword t_n_clust, std::string warping_method, std::string center_method,
    std::string similarity_method, std::string t_optim_method,const rowvec& t_seeds,
    const rowvec& t_warping_opt, const rowvec& t_center_opt, double t_n_out,
    double t_toll, bool t_fence, uword t_iter_max,bool t_show_iter,
    bool t_check_total_similarity, const bool comp_original_center,const  urowvec par_opt)
    : x(t_x), y(t_y),  n_clust(t_n_clust),
      optim_method(t_optim_method), seeds(t_seeds),
      warping_opt(t_warping_opt), n_out(t_n_out),
      toll(t_toll),fence(t_fence), iter_max(t_iter_max),show_iter(t_show_iter),
      check_total_similarity(t_check_total_similarity),com_oc(comp_original_center),parallel_opt(par_opt)
{
    //
    // FActories registration
    //

    // dissimilarity factory
    util::SharedFactory<Dissimilarity> disfac;
    disfac.FactoryRegister<Pearson>("pearson");
    disfac.FactoryRegister<L2>("l2");

    //warping factory
    util::SharedFactory<WarpingFunction> warfac;
    warfac.FactoryRegister<ShiftFunction>("shift");
    warfac.FactoryRegister<DilationFunction>("dilation");
    warfac.FactoryRegister<AffineFunction>("affine");
    warfac.FactoryRegister<NoAlignmentFunction>("noalign");

    //center factory
    util::SharedFactory<CenterMethod> cenfac;
    cenfac.FactoryRegister<Medoid>("medoid");
    cenfac.FactoryRegister<PseudoMedoid>("pseudomedoid");
    cenfac.FactoryRegister<Mean>("mean");

    //Optimizer factory
    util::SharedFactory<OptimizerMethod> optfac;
    optfac.FactoryRegister<Bobyqa>("bobyqa");

    //
    //  check-in
    //

    if(t_show_iter == true)
        {
            Rcpp::Rcout<<"---------------------------------------------"<<endl;
            Rcpp::Rcout<<"Check.in inputs:"<<endl;
        }

    checkIn(x, y, warping_method, center_method, similarity_method,
            optim_method,warfac,disfac,cenfac,optfac,
            warping_opt, t_center_opt, par_opt,t_show_iter);

    optimizer = optfac.instantiate(optim_method);
    dissim =  disfac.instantiate(similarity_method) ;
    warping = warfac.instantiate(warping_method);
    cen = cenfac.instantiate(center_method);
    //parameters center method
    cen->setParameters(t_center_opt);

    n_obs = t_y.n_rows;
    n_camp = t_y.n_cols;
    n_dim = t_y.n_slices;

    if(t_show_iter == true)
        {
            Rcpp::Rcout<<"Loading problem..."<<endl;
            Rcpp::Rcout<<"Number of threads: " << par_opt(0)<<endl;
            Rcpp::Rcout<<"Parallel type: " << par_opt(1)<<endl;
            Rcpp::Rcout<<"Dataset(obs x camp x dim): "<<n_obs<<" x "<<n_camp<<" x "<<n_dim<<endl;
            Rcpp::Rcout<<"Warping Method: "<<warping_method<<endl;
            Rcpp::Rcout<<"Center Method:  "<<center_method<<endl;
            Rcpp::Rcout<<"Dissimilarity Method:  "<<similarity_method<<endl;
            Rcpp::Rcout<<"Optimization Method:  "<< optim_method<<endl;
            Rcpp::Rcout<<"Numbero of cluster: "<<n_clust<<endl;
            Rcpp::Rcout<<"Seeds: ";
            rowvec s=t_seeds+1;
            s.raw_print(Rcpp::Rcout);
        }
}

#ifndef CHECK_IN_HPP
#define CHECK_IN_HPP

#include<RcppArmadillo.h>

bool checkIn(const arma::mat& x,
             const arma::cube& y,
             const std::string warping_method,
             const std::string center_method,
             const std::string similarity_method,
             const std::string optim_method,
             const arma::rowvec& warping_opt,
             const arma::rowvec& center_opt,
             const arma::urowvec par_opt)
{
    Rcpp::Rcout<<"Check-in : ";
    // CONTROLLI SU parallel_option
    if(par_opt(1)==1)
    {
        if(center_method!="medoid")
            Rcpp::stop("Parallelization mode 1 is available only with medoid");
    }
    if(par_opt(1)!=0 & par_opt(1)!=1)
        Rcpp::stop("Parallelization mode could be 0 or 1");

    Rcpp::Rcout<<"Done"<<endl;
    return true;
}
#endif

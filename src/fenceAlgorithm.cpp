#include "fenceAlgorithm.h"
#include "utilityFunctions.h"

void iterativeFence(arma::mat parameters,
                    const arma::uword iter,
                    arma::urowvec& labels,
                    arma::rowvec& index,
                    std::shared_ptr<BaseWarpingFunction>& warping,
                    std::shared_ptr<BaseOptimizerFunction>& optimizer,
                    const arma::cube& templates,
                    const arma::mat& x_reg,
                    const arma::cube& y,
                    const arma::rowvec& x_out,
                    std::shared_ptr<BaseDissimilarityFunction>& dissim,
                    const arma::urowvec& ict,
                    const bool show_iter)
{
    // reading dimension
    arma::uword n_dim = y.n_slices;
    arma::uword n_par = parameters.n_rows;
    arma::urowvec sel_out;
    arma::mat bounds(2,n_par);

    //compute sel_out : observation with outliers
    for(arma::uword i=0; i<n_par; i++)
    {
        double q1 = quantile(parameters.row(i), 0.25);
        double q3 = quantile(parameters.row(i), 0.75);
        bounds(0,i) = q1 - 1.5 * (q3 - q1);
        bounds(1,i) = q3 + 1.5 * (q3 - q1);

        sel_out = join_horiz( sel_out,
                             which_out(parameters.row(i),bounds(0,i),bounds(1,i))
        );
    }

    sel_out = unique(sel_out);
    arma::uword it=1;

    // plot FENCE
    if(show_iter==true)
    {
        Rcpp::Rcout<<"IterativeFence: ";
    }


    while(sel_out.n_cols!=0 && it <10)
    {
        if(show_iter==true)
        {
            Rcpp::Rcout<<" it."<<it<<" ";
        }

        //optimization foe each observation with outliers
        warping->SetParameterBounds(bounds);
        for (unsigned int i = 0;i < sel_out.n_cols;++i)
        {
            unsigned int obs = sel_out[i];
            unsigned int nt = templates.n_rows;
            arma::rowvec index_temp(nt);
            arma::mat parameters_temp(n_par,nt);
            arma::rowvec arg(n_par);
            Rcpp::List interpolationResults = approx(
                x_reg.row(obs),
                y.row(obs),
                "Linear"
            );
            arma::mat y_reg = interpolationResults["values"];

            for (unsigned int t = 0;t < nt;++t)
            {
                arma::mat t_in = templates(arma::span(t), arma::span::all, arma::span::all);

                if (n_dim > 1)
                    t_in = t_in.t();

                WarpingSet wset = warping->SetInputData(x_out, x_out, y_reg, t_in, dissim);

                auto fun = [&] (const arma::rowvec& arg)
                {
                    return warping->GetDissimilarityAfterWarping(wset,arg);
                };

                index_temp(t) = optimizer->Optimize( arg, warping, fun);
                parameters_temp.col(t) = arg;
            }// end iterations on templates

            index(obs) = min( index_temp );
            labels(obs) = ict( index_min(index_temp));
            parameters.col(obs)= parameters_temp.col(index_min(index_temp));
        }//end iterations on sel_out


        //recompute sel out
        sel_out.reset();

        for(arma::uword i=0; i<n_par; i++)
        {
            double q1 = quantile(parameters.row(i),0.25);
            double q3 = quantile(parameters.row(i),0.75);

            bounds(0,i)= q1 - 1.5 * (q3-q1);
            bounds(1,i)= q3 + 1.5 * (q3-q1);

            sel_out= join_horiz( sel_out,
                                 which_out(parameters.row(i),bounds(0,i),bounds(1,i))
            );
        }
        sel_out = unique(sel_out);
        it++;


    }//fine while

    if(show_iter==true)
    {
        Rcpp::Rcout<<" Done."<< std::endl;
    }
    return;
}

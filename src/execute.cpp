#include <RcppArmadillo.h>
#include <Rcpp/Benchmark/Timer.h>

#ifdef _OPENMP
#include <omp.h>
#endif

#include "newcenters.h"
#include "kma_model.h"
#include "optimizer.h"
#include "fence.h"


Rcpp::List KmaModel::execute()
{

    if(show_iter==true)
        cout<<"Start execution."<<endl;

    Rcpp::Timer timer;
    timer.step( "start execution");

    //
    //compute x_out
    //
    if(show_iter==true)  cout<<"Compute x_out: ";
    rowvec x_out = linspace<rowvec>( min(util::lowers(x)), max(util::uppers(x)), n_out);
    if(show_iter==true) cout<<"Done"<<endl;


    //
    //starting template approximated on x_out
    //

    if(show_iter==true) cout<<"Compute initial templates: ";

    field<cube> templates_vec(1,iter_max);
    field<rowvec> x_out_vec(1,iter_max);

    x_out_vec(0)=x_out;
    cube templates(n_clust,n_out,n_dim);

    for(uword i=0; i < n_clust; i++)
        {
            templates(span(i),span::all,span::all) = util::approx( util::abscissa(x,seeds(i)), util::observation(y,seeds(i)), x_out ).t();
        }
    if(show_iter==true) cout<<"Done"<<endl;

    //
    //compute center_origin (to be fixed with new centers)
    //
    if(show_iter==true) cout<<"Compute center_origin and dissimilarity with others : ";

    center original_center;
    if(com_oc==TRUE)
        {
            original_center = cen->computeCenter( x, y,dissim, x_out);
            if(show_iter==true) cout<<"Done"<<endl;
        }
    else
        {
            if(show_iter==true) cout<<"Skipped"<<endl;
        }

    //
    // WHILE equipment
    //
    if(show_iter==true) cout<<"Start while iteration"<<endl;

    // indici di similarità(distanza) ad ogni iterazione
    rowvec index(n_obs);
    index += 1000;
    rowvec index_old(n_obs);
    index_old += 10000;

    // inizializzo  vettore per salare parametri ad ogni iterazione
    uword np = warping->n_pars();
    cube parameters_vec(np,n_obs,iter_max);
    mat x_reg = x;

    // flag for total similarity check
    bool still_in = true;

    // labels del cluster di appartenenza ad ogni iterazione
    urowvec labels(n_obs);
    labels.ones();
    urowvec labels_old(n_obs);

    // indici dei cluster correnti
    urowvec ict = linspace<urowvec>( 0, n_clust-1, n_clust);
    uword iter = 0;


    timer.step("seeds and original center" );

    //if n_clust i want to avoid the check on the labels because they don't change
    uword pn_obs =n_obs;
    if(n_clust==1) pn_obs+=1;

    while(  sum( abs(index-index_old) < toll)<n_obs  && //
            (sum( labels == labels_old  ) != pn_obs) && // non considerà il caso in cui i cluster sn uguali ma con diverse etichette
            (still_in == true) &&
            (iter < iter_max))
        {

            iter++;
            if(show_iter==true) cout<< "Iteration num: "<<iter<<endl;


            index_old = index;
            labels_old = labels;
            mat parameters(np, n_obs);

            if(show_iter==true) cout<<"Set bound: ";
            warping->set_bounds(warping_opt,x_reg);
            if(show_iter==true) cout<<"Done"<<endl;

            //compute best warping parameters and assign new labels
            if(show_iter==true) cout<<iter<<". Compute best warping: "<<endl;

#ifdef _OPENMP
            #pragma omp parallel for num_threads(parallel_opt(0))
#endif
            for(uword obs=0; obs<n_obs; obs++)
                {
                    // inizializzo container warp_temp
                    uword nt = templates.n_rows;
                    rowvec index_temp(nt);
                    mat parameters_temp(np,nt);
                    colvec arg(np);
                    mat y_reg = util::approx( x_reg.row(obs), util::observation(y,obs), x_out );
                    //calcolo warping parameters for each templates
                    for(uword t=0; t<nt; t++ )
                        {
                            mat t_in = templates(span(t),span::all,span::all);
                            if(n_dim >1) t_in =t_in.t();
                            warping_set wset = warping->set_function(x_out, x_out, y_reg, t_in, dissim);

                            //LAMBDA FUNCTION
                            auto fun = [this,&wset] (const colvec& arg)
                            {
                                return this->warping->warp(wset,arg);
                            };

                            index_temp(t) = optimizer->optimize( arg, warping, fun);
                            parameters_temp.col(t) = arg;
                        }

                    //fine iterazioni per ogni tempalte
                    index(obs) = min( index_temp );
                    labels(obs) = ict( index_min(index_temp));
                    parameters.col(obs)= parameters_temp.col(index_min(index_temp));
                    //PRINT
                    #pragma omp critical
                    {
                        if(show_iter==true)
                            {
                                // cout<<iter<<". Observation "<<obs<<" has been assigned to cluster "<<labels(obs)<<endl;
                                // cout<<"distances from other templates are "<<endl;
                                // index_temp.raw_print(std::cout);
                                //cout<<"|";
                            }
                    }
                }// fine iterazioni per ogni osservazione


            if(show_iter==true) cout<<"Done"<<endl;
            timer.step( "warping "+ std::to_string(iter) );


            //update current template list
            ict = unique(labels);


            //PRINT
            if(show_iter==true)
                {
                    //cout<<"Computed labels are:"<<endl;
                    //labels.raw_print();
                    Rcpp::Rcout<<"current cluster vector updated"<<endl;
                    ict.print();
                    std::map<uword,uword> mcl = util::tableC(labels);
                    for(auto it = mcl.cbegin(); it != mcl.cend(); ++it)
                        {
                            Rcpp::Rcout <<"cluster num: "<< it->first << " has " << it->second << " elements;" << endl;
                        }
                    //   cout<<"parameters computed are"<<endl<<parameters<<endl;
                }


            if(show_iter==true) Rcpp::Rcout<<"Fence alghoritm : "<<endl;
            if(fence==TRUE)
                {
                    iterativeFence(parameters, iter, labels, index, warping, optimizer, templates,
                                   x_reg, y, x_out,dissim, ict, show_iter);
                    if(show_iter==true) Rcpp::Rcout<<"Done"<<endl;
                }
            else
                {
                    if(show_iter==true) Rcpp::Rcout<<"Failed"<<endl;
                }


// normalizzaziione
            if(show_iter==true) cout<<"Parameters' normalization :";
            warping->normalize(parameters,ict,labels);
            if(show_iter==true) cout<<"Done"<<endl;


// salvo parametri
            parameters_vec(span::all,span::all,span(iter-1))=parameters;


//update x_reg and x_out
            if(show_iter==true) cout<<"Update x_reg and x_out: ";
            x_reg = warping->apply_warping(x_reg,parameters);
            x_out =  linspace<rowvec>( min( util::lowers(x_reg) ), max( util::uppers(x_reg)), n_out);
            if(show_iter==true) cout<<"Done"<<endl;

            x_out_vec(iter)=x_out;

            timer.step( "fence/norm/update "+ std::to_string(iter) );

            //compute new tempaltes
            if(show_iter==true) cout<<"Compute new templates: "<<endl;

            templates_vec(iter-1) = templates;
            templates.set_size(ict.size(),n_out,n_dim);

            newCenters( x_reg, y, x_out, dissim, cen, parallel_opt,
                         templates, ict, labels, show_iter);

            if(show_iter==true)
                {
                    cout<<"Templates updated"<<endl;
                    cout<<"While condition"<<endl<<"dissim cambiata piu di toll: "<<(sum( abs(index-index_old) < toll)<n_obs )
                        <<endl<<"almeno un etichetta cambiata: "<<(sum( labels == labels_old  ) != n_obs)<<endl;
                    cout<<"Check total similarity: ";
                }

            //check total smilarity
            if(check_total_similarity == true)
                {

                    double tot = sum(index);
                    double tot_old = sum(index_old);
                    // sel la distanza totale aumenta
                    if( tot_old < tot )
                        {
                            still_in= false;
                            templates=templates_vec(iter-1);
                            index = index_old;
                            labels=labels_old;
                            x_out=x_out_vec(iter-1);
                            if(show_iter==true) cout<<"Total similarity didn't increase. ";
                        }
                    if(show_iter==true) cout<<"Done"<<endl;
                }
            timer.step( "newtemplates "+ std::to_string(iter) );

        }//fine while

    parameters_vec.resize(np,n_obs,iter);
    templates_vec(iter) = templates;

    if(show_iter==true) cout<<"End while iterations"<<endl;


    //
    // output
    //
    if(show_iter==true) cout<<"Final warping: ";
    mat final_par = warping->final_warping(parameters_vec,labels,ict);
    if(show_iter==true) cout<<"Done"<<endl;


    field<mat> par_vec(iter);
    for(int k = 0; k < iter; k++)
        {
            par_vec(k) = parameters_vec.slice(k);
        }


    Rcpp::NumericVector out1 = Rcpp::wrap(original_center.x_center);
    out1.attr("dim")=R_NilValue;

    Rcpp::NumericVector out2 = Rcpp::wrap(original_center.dissim_whit_origin);
    out2.attr("dim")=R_NilValue;

    Rcpp::NumericVector out3 = Rcpp::wrap(x_out);
    out3.attr("dim")=R_NilValue;

    Rcpp::NumericVector out4 = Rcpp::wrap(index);
    out4.attr("dim")=R_NilValue;

    Rcpp::NumericVector out7 = Rcpp::wrap(labels+1);
    out7.attr("dim")=R_NilValue;


    field<cube> out8 =templates_vec.cols(0,iter);
    field<rowvec> out9 = x_out_vec.cols(0,iter);


    timer.step( "output ");
    if(show_iter==true) cout<<"Output"<<"---------------------------------------------"<<endl;
    return util::ListBuilder()
           .add("iterations", iter)
           .add("n.clust",n_clust)
           .add("x.center.orig",out1)
           .add("y.center.orig",original_center.y_center)
           .add("similarity.origin",out2)
           .add("x.final", x_reg)
           .add("n.clust.final", ict.size())
           .add("x.centers.final", out3)
           .add("y.centers.final",templates)
           .add("templates_vec",out8)
           .add("x_out_vec",out9)
           .add("labels",out7)
           .add("similarity.final",out4)
           .add("parameters.list", par_vec)
           .add("parameters", final_par)
           .add("timer",timer)
           ;

}

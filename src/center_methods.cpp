
#include <RcppArmadillo.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "center_methods.h"
#include "low.h"
#include "utilities.h"

//
//  COMPUTE CENTERS
//


center Medoid::computeCenter(const mat& x,const cube& y, std::shared_ptr<Dissimilarity>& dissim, const rowvec& x_out)
{
    center out;

    out.x_center = x_out;
    uword n_out = x_out.size();
    rowvec x_com = linspace<rowvec>(max(util::lowers(x)),min(util::uppers(x)), n_out);
    uword n_obs = y.n_rows;
    uword n_dim = y.n_slices;

    mat y_fin(n_dim,n_out);


    mat D(n_obs,n_obs);
    D.zeros();

    for(uword i=0; i<n_obs; i++)
        {
            for(uword j=(i+1); j<n_obs; j++)
                {

                    mat obs_i = y(span(i),span::all,span::all);
                    mat obs_j = y(span(j),span::all,span::all);
                    if(n_dim >1)
                        {
                            obs_i = obs_i.t();
                            obs_j = obs_j.t();
                        }

                    D(i,j) = dissim->compute( x_com, x_com,
                                              util::approx( x.row(i), obs_i, x_com),
                                              util::approx( x.row(j), obs_j, x_com));
                    D(j,i) = D(i,j);
                }
        }

    colvec dis = sum(D,1);
    uword m =  index_min(dis);
    mat obs_m = y(span(m),span::all,span::all);

    if(n_dim >1) obs_m = obs_m.t();

    out.y_center =  util::approx( x.row(m), obs_m, out.x_center);

    //
    // compute dissimilarity whit others
    //
    out.dissim_whit_origin = D.row(m);

    return out;
}

center Mean::computeCenter(const mat& x,const cube& y, std::shared_ptr<Dissimilarity>& dissim, const rowvec& x_out)
{

    center out;

    uword n_out = x_out.size();
    uword n_dim = y.n_slices;
    uword n_obs = y.n_rows;
    uword n_cam = y.n_cols;
    mat y_fin(n_dim,n_out);
    for(uword l=0; l<n_dim; l++)
        {

            vector<double> y_in, x_in;

            for(uword i=0; i<n_obs; i++)
                {
                    for(uword j=0; j<n_cam; j++)
                        {

                            if ( is_finite(x(i,j)) &&is_finite(y(i,j,l)))
                                {
                                    x_in.push_back( x(i,j) );
                                    y_in.push_back( y(i,j,l) );
                                }
                        }
                }// fine ciclu su osservaizoni

            std::vector<std::pair<double,double>> zipped;
            util::zip(x_in, y_in, zipped);

            // Sort the vector of pairs
            std::sort(std::begin(zipped), std::end(zipped),
                      [&](const pair<double,double>& a, const pair<double,double>& b)
            {
                return a.first < b.first;
            });

            util::unzip(zipped,x_in,y_in);
            vector<double> ys(y_in.size());
            lowess(x_in, y_in, span, d, 2, ys);

            rowvec x1(x_in),y1(ys);
            // approssimo su x_out
            y_fin.row(l) = util::approx(x1,y1,x_out);

        }
    out.y_center= y_fin;
    // assegno x_center
    out.x_center=x_out;

    // compute dissimilarity whit others
    rowvec dso(n_obs);
    for(uword i=0; i<n_obs; i++)
        {
            dso(i)= dissim->compute(x_out,x.row(i),out.y_center,util::observation(y,i));
        }
    out.dissim_whit_origin=dso;

    return out;
}

center PseudoMedoid::computeCenter(const mat& x,const cube& y, std::shared_ptr<Dissimilarity>& dissim, const rowvec& x_out)
{

    center out;

    out.x_center = x_out;
    uword n_out = x_out.size();

    rowvec x_com = linspace<rowvec>(max(util::lowers(x)),min(util::uppers(x)), n_out);

    uword n_obs = y.n_rows;
    uword n_dim = y.n_slices;
    mat y_fin(n_dim,n_out);

    for(uword k=0; k<n_dim; k++)
        {

            mat D(n_obs,n_obs); // D.zeros();
            mat dim = y.slice(k);

            for(uword i=0; i<n_obs; i++)
                {
                    for(uword j=i; j<n_obs; j++)
                        {

                            D(i,j) = dissim->compute( x_com, x_com,
                                                      util::approx( x.row(i), dim.row(i), x_com),
                                                      util::approx( x.row(j), dim.row(j), x_com));
                            D(j,i) = D(i,j);
                        }
                }

            colvec dis = sum(D,1);
            uword m =  index_min(dis);
            y_fin.row(k) =  util::approx( x.row(m), dim.row(m), out.x_center);
        }


    out.y_center=util::norm_ex(y_fin);

    //
    // compute dissimilarity whit others
    //
    rowvec dso(n_obs);
    for(uword i=0; i<n_obs; i++)
        {
            dso(i)= dissim->compute(out.x_center, x.row(i), out.y_center, util::observation(y,i) );
        }
    out.dissim_whit_origin=dso;
    return out;
}

//
// parallel center
//
center Medoid::computeParallelCenter(const mat& x, const cube& y, std::shared_ptr<Dissimilarity>& dissim, const rowvec& x_out,uword n_th)
{
    center out;

    out.x_center = x_out;
    uword n_out = x_out.size();
    rowvec x_com = linspace<rowvec>(max(util::lowers(x)),min(util::uppers(x)), n_out);
    uword n_obs = y.n_rows;
    uword n_dim = y.n_slices;

    mat y_fin(n_dim,n_out);

    field<rowvec> fD(n_obs);

    for(uword i=0; i<n_obs;i++){
      fD(i).zeros(n_obs);
    }

#ifdef _OPENMP
    #pragma omp parallel for num_threads(n_th)
#endif
    for(uword k=1; k<= n_obs*(n_obs-1)/2 ; k++)
        {
            double kd =k;
            double i = floor( (1+sqrt(8*kd - 7)) / 2 );
            double j = kd-(i-1)*i/2-1;


            mat obs_i = y(span(i),span::all,span::all);
            mat obs_j = y(span(j),span::all,span::all);

            if(n_dim >1)
                {
                    obs_i = obs_i.t();
                    obs_j = obs_j.t();
                }


            // D(i,j)
            fD(i)(j)= dissim->compute( x_com, x_com,
                                      util::approx( x.row(i), obs_i, x_com),
                                      util::approx( x.row(j), obs_j, x_com));
            fD(j)(i) = fD(i)(j);
        }

    //colvec dis = sum(D,1);
    colvec dis(n_obs);
    for(uword i=0; i<n_obs;i++)
      dis(i) = sum(fD(i));

    uword m =  index_min(dis);
    mat obs_m = y(span(m),span::all,span::all);

    if(n_dim >1) obs_m = obs_m.t();

    out.y_center =  util::approx( x.row(m), obs_m, out.x_center);

    //
    // compute dissimilarity whit others
    //
    out.dissim_whit_origin = fD(m);

    return out;
}


//
// set paramenters
//
void Mean::setParameters(const rowvec& par)
{
    span = par(0);
    d = par(1);
}

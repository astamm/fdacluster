#include "warping.h"

//
// NORMALIZE
//

void AffineFunction::normalize(mat& par,const urowvec& ict,const urowvec& labels)
{

    for(size_t i=0; i<ict.size(); i++)
    {

        urowvec sel =find(labels==ict(i)).t();
        //calcolo medie cluster
        colvec par_mean= mean(par.cols(sel),1);

        // aggiorno shift e dilation
        for(size_t j=0; j<sel.size(); j++)
        {
            //normalized dilation
            par(0,sel(j))=  par(0,sel(j))/par_mean(0);
            //normalized shift
            par(1,sel(j))= -par(0,sel(j)) * par_mean(1)/par_mean(0) + par(1,sel(j));
        }
    }// for su ict
}
void DilationFunction::normalize(mat& par,const urowvec& ict,const urowvec& labels)
{

    for(size_t i=0; i<ict.size(); i++)
    {
        urowvec sel =find(labels==ict(i)).t();

        //calcolo medie cluster
        colvec par_mean= mean(par.cols(sel),1);

        // aggiorno shift e dilation
        for(size_t j=0; j<sel.size(); j++)
        {
            par(0,sel(j))=  par(0,sel(j))/par_mean(0);
        }// for su sel

    }// for su ict
}
void ShiftFunction::normalize(mat& par,const urowvec& ict,const urowvec& labels)
{

    for(size_t i=0; i<ict.size(); i++)
    {
        urowvec sel =find(labels==ict(i)).t();

        //calcolo medie cluster
        colvec par_mean= mean(par.cols(sel),1);

        // aggiorno shift e dilation
        for(size_t j=0; j<sel.size(); j++)
        {

            par(0,sel(j))=  par(0,sel(j)) -par_mean(0);

        }// for su sel
    }// for su ict
}

//
// Final Warping
//
mat AffineFunction::final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict)
{
    uword np = parameters_vec.n_rows;
    uword n_obs = parameters_vec.n_cols;
    uword n_iter = parameters_vec.n_slices;
    mat out(np,n_obs);
    out.row(0).ones();
    for(size_t i=0; i<n_iter; i++)
    {
        rowvec a = parameters_vec(span(0),span::all,span(i));
        rowvec b = parameters_vec(span(1),span::all,span(i));
        out.row(0) = out.row(0) % a;
        out.row(1) = out.row(1) % a + b;
    }
    for(uword k = 0; k < ict.size(); k++)
    {
        urowvec sel = find(labels == ict(k)).t();
        //compute means
        colvec m = sum(out.cols(sel),1)/sel.size();
        for(size_t i =0; i< sel.size(); i++)
        {
            out(0,sel(i)) = out(0,sel(i))/m(0);
            out(1,sel(i)) = (out(1,sel(i)) - m(1))/m(0);
        }
    }
    return out;
}
mat ShiftFunction::final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict)
{
    uword n_obs = parameters_vec.n_cols;
    uword n_iter = parameters_vec.n_slices;
    mat out(1,n_obs);
    out.row(0).zeros();
    for(size_t i=0; i<n_iter; i++)
    {
        rowvec a = parameters_vec(span(0),span::all,span(i));
        out.row(0) = out.row(0) + a;
    }
    for(uword k = 0; k < ict.size(); k++)
    {
        urowvec sel = find(labels == ict(k)).t();
        colvec m = mean(out.cols(sel),1);
        for(uword i =0; i< sel.size(); i++)
        {
            out(0,sel(i)) = (out(0,sel(i)) - m(0));
        }
    }
    return out;
}
mat DilationFunction::final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict)
{
    uword n_obs = parameters_vec.n_cols;
    uword n_iter = parameters_vec.n_slices;
    mat out(1,n_obs);
    out.row(0).ones();
    for(size_t i=0; i<n_iter; i++)
    {
        rowvec a = parameters_vec(span(0),span::all,span(i));
        out.row(0) = out.row(0) % a ;
    }
    for(uword k = 0; k < ict.size(); k++)
    {
        urowvec sel = find(labels == ict(k)).t();
        //compute means
        colvec m = mean(out.cols(sel),1);
        for(size_t i =0; i< sel.size(); i++)
        {
            out(0,sel(i)) = out(0,sel(i))/m(0);
        }
    }
    return out;
}
mat NoAlignmentFunction::final_warping(const cube& parameters_vec, const urowvec& labels, const urowvec& ict)
{
    mat out(0,labels.n_cols);
    return out;

}


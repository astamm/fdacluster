#include <RcppArmadillo.h>

#include "utilities.h"

using namespace arma;

//
// tableC
//

std::map<arma::uword, arma::uword> util::tableC(arma::urowvec x)
{

    std::map<arma::uword,arma::uword> counts;
    int n = x.size();

    for (int i = 0; i < n; i++)
        {
            counts[x[i]]=counts[x[i]]+1;
        }

    return counts;
}








//
//  which_out
//

urowvec util::which_out(rowvec vec, double lo,double up)
{

    urowvec out(vec.n_cols);
    uword c=0;
    for(int i=0; i<vec.size(); i++)
        {
            if(vec(i) < lo || vec(i) > up)
                {
                    out(c)=i;
                    c++;
                }
        }
    out.resize(c);
    return out;
}



//
//  quantile
//


double util::quantile(const arma::rowvec vec, const double per)
{

    if(per<=0 || per>=1)
        perror("quantile:per has to be between 0 and 1");

    std::vector<double> svec = conv_to<std::vector<double>>::from(vec);
    uword q = round(vec.n_cols * per) ;

    std::nth_element(svec.begin(),svec.begin() + q, svec.end());

    return *(svec.begin()+q);

}



//
//  norm_estremi
//

mat util::norm_ex(const mat& y)
{
    mat out(y.n_rows,y.n_cols);
    out.fill(datum::nan);

    for(size_t i=0; i< y.n_cols; i++)
        {
            if(y.col(i).is_finite())
                {
                    out.col(i)=y.col(i);
                }
        }

    return out;
}



//
// Upper and Lower element of each row not NA
//
//
rowvec util::uppers(const mat& x)
{
    rowvec out(x.n_rows);

    for(size_t i=0; i< x.n_rows; i++)
        {
            out[i]=max(x.row(i));
        }
    return out;
}

rowvec util::lowers(const mat& x)
{
    rowvec out(x.n_rows);

    for(size_t i=0; i< x.n_rows; i++)
        {
            out[i]=min(x.row(i));
        }
    return out;
}

//
// observation
//
const mat util::observation(const cube& y,uword i)
{

    if(y.n_slices > 1)
        {
            mat out =y( span(i,i), span::all, span::all );
            return out.t();
        }
    return y( span(i,i), span::all, span::all );
}


const cube util::observations(const cube& y, urowvec ind)
{
    cube out(ind.size(),y.n_cols,y.n_slices);

    for(uword i =0; i< ind.size(); i++)
        {
            out( span(i,i), span::all, span::all) = y( span(ind(i),ind(i)), span::all, span::all );
        }
    return out;
}

//
// ABSCISSA
//

//
// abscissa
//
const rowvec util::abscissa(const mat& x,uword i)
{
    return x.row(i);
}

//
// APPROX
//
mat util::approx(const rowvec& x,
                 const mat& y,
                 const rowvec& xx)
{

    uword dim = y.n_rows;
    uword n_xx  = xx.n_cols;

    mat yy(dim, n_xx);
    yy.fill(datum::nan);

    //check sul vettore in ingresso non necessaria ma utile nello sviluppo
    if(xx(n_xx-1)-xx(0) <= 0)
        {
            cout<< " diff :"<<xx(n_xx-1)-xx(0) <<endl;
            cout<<" vettore in ingresso decrescente"<<endl;
            return yy;
        }

    uword i=0;

    while(xx(i) <= x(0))
        {
            if(xx(i) == x(0))
                {
                    // for(uword k = 0; k < dim; k++)
                    //     {
                    //         yy(k,i)=y(k,0);
                    //     }
                    yy.col(i)=y.col(0);
                }
            i++;
        }

    for(uword j = 1; j < x.size(); j++)
        {

            //  il <= non funziona bene con double quindi per l'uguaglianza controll la diff< eps
            while( ( i < xx.size() ) && (xx(i) < x(j) || (std::abs(xx(i) - x(j)) < 0.000000001) ) )
                {
                    // for(int k = 0; k < dim; k++)
                    //     {
                    //         yy(k,i)= ( y(k,j) - y(k,j-1) ) * ( xx(i) - x(j-1)) / (x(j)-x(j-1)) + y(k,j-1);
                    //     }
                    yy.col(i)= y.col(j-1) +( y.col(j) - y.col(j-1) ) * ( xx(i) - x(j-1)) / (x(j)-x(j-1)) ;
                    i++;
                }

        }

    return yy;

}


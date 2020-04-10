// Copyright (C) 2017 Alessandro Zito (zito.ales@gmail.com)
//
// This file is part of Fdakmapp.
//
// Fdakmapp is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// Fdakmapp is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
// GNU General Public License for more details.
//
//   You should have received a copy of the GNU General Public License
//   along with Fdakmapp.  If not, see <http://www.gnu.org/licenses/>.

#include "lowess.h"

#include <algorithm>
#include <cmath>
#include <fstream>

void lowess(const std::vector<double> &x,
            const std::vector<double> &y,
            double f,
            double delta,
            long nsteps,
            std::vector<double> &ys)
{
    std::vector<double> rw;
    std::vector<double> res;
    lowess(x, y, f, nsteps, delta, ys, rw, res);
}

void lowess(const std::vector<double> &x,
            const std::vector<double> &y,
            double f,
            long nsteps,
            double delta,
            std::vector<double> &ys,
            std::vector<double> &rw,
            std::vector<double> &res)
{
    long n=(long)x.size();
    bool ok=false;
    long nleft,nright, i, j, iter, last, m1, m2, ns;
    double cut, cmad, r, d1, d2, c1, c9, alpha, denom;
    if((n==0)||((long)y.size()!=n)) return;
    ys.resize(n);
    rw.resize(n);
    res.resize(n);
    if(n==1)
    {
        ys[0]=y[0];
        return;
    }
    // ns - at least 2, at most n
    ns = std::max(std::min((long)(f*n),n),(long)2);
    for(iter=0; iter<nsteps+1; iter++)
    {
        // robustness iterations
        nleft = 0;
        nright = ns-1;
        // index of last estimated point
        last = -1;
        // index of current point
        i=0;
        do
        {
            while(nright<n-1)
            {
                // move <nleft,nright> right, while radius decreases
                d1 = x[i]-x[nleft];
                d2 = x[nright+1] - x[i];
                if(d1<=d2)break;
                nleft++;
                nright++;
            }
            // fit value at x[i]
            lowest(x,y,x[i],ys[i],nleft,nright,res,iter>0,rw,ok);
            if(!ok) ys[i]=y[i];
            if(last<i-1)
            {
                // interpolate skipped points
                if(last<0)
                {
                    //warning("Lowess: out of range.\n");
                }
                denom = x[i] - x[last];
                for(j=last+1; j<i; j++)
                {
                    alpha = (x[j]-x[last])/denom;
                    ys[j] = alpha * ys[i] + (1.0-alpha)*ys[last];
                }
            }
            last = i;
            cut = x[last]+delta;
            for(i=last+1; i<n; i++)
            {
                if(x[i]>cut)break;
                if(x[i]==x[last])
                {
                    ys[i]=ys[last];
                    last=i;
                }
            }
            i=std::max(last+1,i-1);
        }
        while(last<n-1);
        for(i=0; i<n; i++)
            res[i] = y[i]-ys[i];
        if(iter==nsteps)break ;
        for(i=0; i<n; i++)
            rw[i]=abs(res[i]);
        sort(rw.begin(),rw.end());
        m1 = n/2+1;
        m2 = n-m1;
        m1 --;
        cmad = 3.0 *(rw[m1]+rw[m2]);
        c9 = .999*cmad;
        c1 = .001*cmad;
        for(i=0; i<n; i++)
        {
            r = abs(res[i]);
            if(r<=c1) rw[i]=1;
            else if(r>c9) rw[i]=0;
            else rw[i] = (1.0-(r/cmad)*(r/cmad))*(1.0-(r/cmad)*(r/cmad));
        }
    }
}

void lowest(const std::vector<double> &x,
            const std::vector<double> &y,
            double xs,
            double &ys,
            long nleft,
            long nright,
            std::vector<double> &w,
            bool userw,
            std::vector<double> &rw,
            bool &ok)
{
    long n = (long)x.size();
    long nrt, j;
    double a, b, c, h, r, h1, h9, range;
    range = x[n-1]-x[0];
    h = std::max(xs-x[nleft],x[nright]-xs);
    h9 = 0.999*h;
    h1 = 0.001*h;
    // sum of weights
    a = 0;
    for(j=nleft; j<n; j++)
    {
        // compute weights (pick up all ties on right)
        w[j]=0.;
        r = abs(x[j]-xs);
        if(r<=h9)
        {
            // small enough for non-zero weight
            if(r>h1) w[j] = (1.0-(r/h)*(r/h)*(r/h))*(1.0-(r/h)*(r/h)*(r/h))*(1.0-(r/h)*(r/h)*(r/h));
            else w[j] = 1.;
            if(userw) w[j] *= rw[j];
            a += w[j];
        }
        else if(x[j]>xs) break;  // get out at first zero wt on right
    }
    nrt = j-1;
    // rightmost pt (may be greater than nright because of ties)
    if(a<=0.) ok = false;
    else
    {
        // weighted least squares
        ok = true;
        // normalize weights
        for(j=nleft; j<=nrt; j++)
            w[j] /= a;
        if(h>0.)
        {
            // use linear fit
            a = 0.;
            for(j=nleft; j<=nrt; j++)
                a += w[j]*x[j]; // weighted centre of values
            b = xs-a;
            c = 0;
            for(j=nleft; j<=nrt; j++)
                c += w[j]*(x[j]-a)*(x[j]-a);
            if(sqrt(c)>0.001*range)
            {
                // points are spread enough to compute slope
                b /= c;
                for(j=nleft; j<=nrt; j++)
                    w[j] *= (1.0+b*(x[j]-a));
            }
        }
        ys = 0;
        for(j=nleft; j<=nrt; j++)
            ys += w[j]*y[j];
    }
}

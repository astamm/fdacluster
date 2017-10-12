
#ifndef LOW_H
#define LOW_H

#include<vector>

using namespace std;
/// Lowess algorithm
/**
 * Lowess algorithm for the omputation of new centers in k-mean algorithm.
 */
void lowess(const std::vector<double> &x, const std::vector<double> &y, double f, long nsteps, double delta, std::vector<double> &ys, std::vector<double> &rw, std::vector<double> &res);

void lowess(const std::vector<double> &x, const std::vector<double> &y, double f,double delta, long nsteps, std::vector<double> &ys);

void lowest(const std::vector<double> &x, const std::vector<double> &y, double xs, double &ys, long nleft, long nright, std::vector<double> &w,bool userw,  std::vector<double> &rw, bool &ok);

#endif

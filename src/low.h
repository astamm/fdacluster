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

#ifndef LOW_H
#define LOW_H

#include <vector>

using namespace std;
/// Lowess algorithm
/**
 * Lowess algorithm for the omputation of new centers in k-mean algorithm.
 */
void lowess(const std::vector<double> &x, const std::vector<double> &y, double f, long nsteps, double delta, std::vector<double> &ys, std::vector<double> &rw, std::vector<double> &res);

void lowess(const std::vector<double> &x, const std::vector<double> &y, double f,double delta, long nsteps, std::vector<double> &ys);

void lowest(const std::vector<double> &x, const std::vector<double> &y, double xs, double &ys, long nleft, long nright, std::vector<double> &w,bool userw,  std::vector<double> &rw, bool &ok);

#endif

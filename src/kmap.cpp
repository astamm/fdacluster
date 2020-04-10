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

#include "kmap.h"

#include "dissimilarity.h"
#include "warping.h"
#include "utilities.h"
#include "kma_model.h"

Rcpp::List kmap(const arma::mat &x,
                const arma::cube &y,
                const arma::urowvec &seeds,
                const unsigned int &n_clust,
                const unsigned int &maximum_number_of_iterations,
                const unsigned int &number_of_threads,
                const unsigned int &parallel_method,
                const double &shift_upper_bound,
                const double &dilation_upper_bound,
                const double &tolerance,
                const double &lowess_span_value,
                const bool &use_fence,
                const bool &check_total_similarity,
                const bool &use_verbose,
                const bool &compute_original_centers,
                const std::string &interpolation_method,
                const std::string &warping_method,
                const std::string &center_method,
                const std::string &dissimilarity_method,
                const std::string &optimizer_method)
{
    KmaModel model;

    model.SetInputData(x, y);

    model.SetSeedVector(seeds);

    model.SetNumberOfClusters(n_clust);
    model.SetMaximumNumberOfIterations(maximum_number_of_iterations);
    model.SetNumberOfThreads(number_of_threads);
    model.SetParallelMethod(parallel_method);

    model.SetShiftUpperBound(shift_upper_bound);
    model.SetDilationUpperBound(dilation_upper_bound);
    model.SetTolerance(tolerance);

    model.SetUseFence(use_fence);
    model.SetCheckTotalSimilarity(check_total_similarity);
    model.SetUseVerbose(use_verbose);
    model.SetComputeOriginalCenters(compute_original_centers);

    model.SetInterpolationMethod(interpolation_method);
    model.SetWarpingMethod(warping_method);
    model.SetCenterMethod(center_method, lowess_span_value);
    model.SetDissimilarityMethod(dissimilarity_method);
    model.SetOptimizerMethod(optimizer_method);

    if (use_verbose)
        model.Print(warping_method, center_method, dissimilarity_method, optimizer_method);

    return model.FitModel();
}

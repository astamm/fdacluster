#include "kmap.h"
#include "kmaModelClass.h"

Rcpp::List kmap(const arma::mat &x,
                const arma::cube &y,
                const arma::urowvec &seeds,
                const arma::rowvec &warping_options,
                const unsigned int &n_clust,
                const unsigned int &maximum_number_of_iterations,
                const unsigned int &number_of_threads,
                const unsigned int &parallel_method,
                const double &distance_relative_tolerance,
                const bool &use_fence,
                const bool &check_total_dissimilarity,
                const bool &use_verbose,
                const bool &compute_overall_center,
                const std::string &warping_method,
                const std::string &center_method,
                const std::string &dissimilarity_method,
                const std::string &optimizer_method)
{
    KmaModel model;

    model.SetInputData(x, y);

    model.SetSeedVector(seeds);
    model.SetWarpingOptions(warping_options);

    model.SetNumberOfClusters(n_clust);
    model.SetMaximumNumberOfIterations(maximum_number_of_iterations);
    model.SetNumberOfThreads(number_of_threads);
    model.SetParallelMethod(parallel_method);

    model.SetDistanceRelativeTolerance(distance_relative_tolerance);

    model.SetUseFence(use_fence);
    model.SetCheckTotalDissimilarity(check_total_dissimilarity);
    model.SetUseVerbose(use_verbose);
    model.SetComputeOverallCenter(compute_overall_center);

    model.SetWarpingMethod(warping_method);
    model.SetCenterMethod(center_method);
    model.SetDissimilarityMethod(dissimilarity_method);
    model.SetOptimizerMethod(optimizer_method);

    if (use_verbose)
        model.Print(warping_method, center_method, dissimilarity_method, optimizer_method);

    return model.FitModel();
}

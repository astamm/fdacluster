#' K-mean alignment and variants for functional data
#'
#' @param x A matrix of size nObs x nPts storing the evaluation grid of each
#'   observation.
#' @param y An 3D array of size nObs x nDim x nPts storing the observation
#'   values.
#' @param seeds A vector of integers of size \code{n_clust} specifying the
#'   indices of the initial templates. Defaults to \code{NULL}, which boils down
#'   to randomly sampled indices.
#' @param warping_options A numeric vector supplied as a helper to the chosen
#'   \code{warping_method} to decide on warping parameter bounds.
#' @param n_clust An integer specifying the number of clusters (default: 1).
#' @param maximum_number_of_iterations An integer specifying the maximum number
#'   of iterations before the algorithm stops (default: 100).
#' @param number_of_threads An integer specifying the number of threads used for
#'   parallelization (default: 1).
#' @param parallel_method An integer value specifying the type of desired
#'   parallelization for template computation, If 0 (default), templates are
#'   computed in parallel. If 1, parallelization occurs within a single template
#'   computation (only for the medoid method as of now).
#' @param distance_relative_tolerance A number specifying a relative tolerance
#'   on the distance update between two iterations. If all observations have not
#'   sufficiently improved in that sense, the algorithm stops. Defaults to 1e-3.
#' @param use_fence A boolean specifying whether the fence algorithm should be
#'   used to robustify the algorithm against outliers (default: \code{FALSE}).
#' @param check_total_dissimilarity A boolean specifying whether an additional
#'   stopping criterion based on improvement of the total dissimilarity should
#'   be used (default: \code{TRUE}).
#' @param use_verbose A boolean specifying whether the algorithm should output
#'   details of the steps to the console (default: \code{TRUE}).
#' @param compute_overall_center A boolean specifying whether the overall center
#'   should be also computed (default: \code{FALSE}).
#' @param warping_method A string specifying the warping method. Choices are
#'   \code{"none"}, \code{"shift"}, \code{"dilation"} and \code{"affine"}
#'   (default).
#' @param center_method A string specifying the center method. Choices are
#'   \code{"medoid"} and \code{"mean"} (default).
#' @param dissimilarity_method A string specifying the dissimilarity method.
#'   Choices are \code{"pearson"} and \code{"l2"} (default).
#' @param optimizer_method A string specifying the optimizer method. The only
#'   choice for now is \code{"bobyqa"}.
#'
#' @return The function output is a \code{kmap} object, which is a list with the
#'   following elements:
#'   \item{x}{As input.}
#'   \item{y}{As input.}
#'   \item{seeds}{Indices used in the algorithm.}
#'   \item{iterations}{Number of iterations before the KMA algorithm stops.}
#'   \item{n_clust}{As input.}
#'   \item{overall_center_grid}{Overall center grid if
#'   \code{compute_overall_center} is set.}
#'   \item{overall_center_values}{Overall center values if
#'   \code{compute_overall_center} is set.}
#'   \item{distances_to_overall_center}{Distances of each observation to the
#'   overall center if \code{compute_overall_center} is set.}
#'   \item{x_final}{Aligned observation grids.}
#'   \item{n_clust_final}{Final number of clusters. Note that
#'   \code{n_clust_final} may differ from initial number of clusters
#'   \code{n_clust} if some clusters are empty.}
#'   \item{x_centers_final}{Final center grids.}
#'   \item{y_centers_final}{Final center values.}
#'   \item{template_grids}{List of template grids at each iteration.}
#'   \item{template_values}{List of template values at each iteration.}
#'   \item{labels}{Cluster memberships.}
#'   \item{final_dissimilarity}{Distances of each observation to the center of
#'   its assigned cluster.}
#'   \item{parameters_list}{List of estimated warping parameters at each
#'   iteration.}
#'   \item{parameters}{Final estimated warping parameters.}
#'   \item{timer}{Execution time step by step.}
#'   \item{warping_method}{As input.}
#'   \item{dissimilarity_method}{As input.}
#'   \item{center_method}{As input.}
#'   \item{optimizer_method}{As input.}
#'
#' @export
#'
#' @examples
#' res <- kma(
#'   simulated30$x,
#'   simulated30$y,
#'   seeds = c(1, 21),
#'   n_clust = 2,
#'   center_method = "medoid",
#'   warping_method = "affine",
#'   dissimilarity_method = "pearson"
#' )
kma <- function(x, y,
                seeds = NULL,
                warping_options = c(0.15, 0.15),
                n_clust = 1,
                maximum_number_of_iterations = 100,
                number_of_threads = 1,
                parallel_method = 0,
                distance_relative_tolerance = 0.001,
                use_fence = FALSE,
                check_total_dissimilarity = TRUE,
                use_verbose = TRUE,
                compute_overall_center = FALSE,
                warping_method = "affine",
                center_method = "mean",
                dissimilarity_method = "l2",
                optimizer_method = "bobyqa") {
  if (anyNA(x))
    stop("The input parameter x should not contain non-finite values.")

  if (anyNA(y))
    stop("The input parameter y should not contain non-finite values.")

  # Handle one-dimensional data
  if (length(dim(y)) == 2) {
    y <- array(y, c(dim(y)[1], 1, dim(y)[2]))
  }

  # Handle vector grid
  if (is.vector(x)) {
    x <- matrix(x, dim(y)[1], dim(y)[3], byrow = TRUE)
  }

  # Handle seeds
  if (is.null(seeds)) {
    seeds <- sample(0:(dim(y)[1] - 1), n_clust)
  } else {
    seeds <- seeds - 1
  }

  out <- kmap(
    x,
    y,
    seeds,
    warping_options,
    n_clust,
    maximum_number_of_iterations,
    number_of_threads,
    parallel_method,
    distance_relative_tolerance,
    use_fence,
    check_total_dissimilarity,
    use_verbose,
    compute_overall_center,
    warping_method,
    center_method,
    dissimilarity_method,
    optimizer_method
  )

  ## gestione timer  ################################################
  time <- diff(round(out$timer / 1000000000, 4))
  t <- data.frame(0, 0, 0, 0, 0, 0)
  names(t) <-
    c("start", "warping", "fece/norm", "templates", "output", "total")
  rownames(t) <- c("sec")
  t[1] <- time[1]
  for (i in 0:(out$iterations - 1)) {
    t[2] <- t[2] + time[2 + (i * 3)]
    t[3] <- t[3] + time[3 + (i * 3)]
    t[4] <- t[4] + time[4 + (i * 3)]
  }
  t[5] <- time[out$iterations * 3 + 2]
  t[6] <- out$timer[length(out$timer)] / 1000000000
  out$timer <- t
  #######################################################################

  out <- c(
    out,
    warping_method = list(warping_method),
    dissimilarity_method = list(dissimilarity_method),
    center_method = list(center_method),
    optimizer_method = list(optimizer_method)
  )

  class(out) <- "kma"

  out
}

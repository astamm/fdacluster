#' K-Means Clustering for Functional Data with Amplitude and Phase Separation
#'
#' This function provides implementations of the k-means clustering algorithm
#' for functional data, with possible joint amplitude and phase separation. A
#' number of warping class are implemented to achieve this separation.
#'
#' @param x A numeric matrix of shape \eqn{N \times M} specifying the grids of
#'   size \eqn{M} on which each of the \eqn{N} curves have been observed.
#' @param y A numeric array of shape \eqn{N \times L \times M} specifying the
#'   \eqn{N}-sample of \eqn{L}-dimensional curves observed on grids of size
#'   \eqn{M}.
#' @param n_clusters An integer value specifying the number of clusters.
#'   Defaults to `1L`.
#' @param warping_class A string specifying the warping class Choices are
#'   `"affine"`, `"dilation"`, `"none"`, `"shift"` or `"srsf"`. Defaults to
#'   `"affine"`. The SRSF class is the only class which is boundary-preserving.
#' @param seeds An integer vector of length `n_clusters` specifying the indices
#'   of the initial templates. Defaults to `NULL`, which boils down to randomly
#'   sampled indices.
#' @param centroid_type A string specifying the type of centroid to compute.
#'   Choices are `"mean"` or `"medoid"`. Defaults to `"mean"`.
#' @param maximum_number_of_iterations An integer specifying the maximum number
#'   of iterations before the algorithm stops if no other convergence criterion
#'   was met. Defaults to `100L`.
#' @param use_verbose A boolean specifying whether the algorithm should output
#'   details of the steps to the console. Defaults to `TRUE`.
#' @param metric A string specifying the metric used to compare curves.
#'   Choices are `"l2"` or `"pearson"`. Defaults to `"l2"`. This is used only
#'   when `warping_class != "srsf"`.
#' @param warping_options A numeric vector supplied as a helper to the chosen
#'   `warping_class` to decide on warping parameter bounds. This is used only
#'   when `warping_class != "srsf"`.
#' @param number_of_threads An integer value specifying the number of threads
#'   used for parallelization. Defaults to `1L`. This is used only when
#'   `warping_class != "srsf"`.
#' @param parallel_method An integer value specifying the type of desired
#'   parallelization for template computation, If `0L`, templates are computed
#'   in parallel. If `1L`, parallelization occurs within a single template
#'   computation (only for the medoid method as of now). Defaults to `0L`. This
#'   is used only when `warping_class != "srsf"`.
#' @param distance_relative_tolerance A numeric value specifying a relative
#'   tolerance on the distance update between two iterations. If all
#'   observations have not sufficiently improved in that sense, the algorithm
#'   stops. Defaults to `1e-3`. This is used only when `warping_class !=
#'   "srsf"`.
#' @param use_fence A boolean specifying whether the fence algorithm should be
#'   used to robustify the algorithm against outliers. Defaults to `FALSE`. This
#'   is used only when `warping_class != "srsf"`.
#' @param check_total_dissimilarity A boolean specifying whether an additional
#'   stopping criterion based on improvement of the total dissimilarity should
#'   be used. Defaults to `TRUE`. This is used only when `warping_class !=
#'   "srsf"`.
#' @param compute_overall_center A boolean specifying whether the overall center
#'   should be also computed. Defaults to `FALSE`. This is used only when
#'   `warping_class != "srsf"`.
#'
#' @return An object of class [`caps`].
#'
#' @export
#' @examples
#' res <- fdakmeans(
#'   simulated30$x,
#'   simulated30$y,
#'   seeds = c(1, 21),
#'   n_clusters = 2,
#'   centroid_type = "medoid",
#'   warping_class = "affine",
#'   metric = "pearson"
#' )
fdakmeans <- function(x, y,
                      n_clusters = 1L,
                      warping_class = c("affine", "dilation", "none", "shift", "srsf"),
                      seeds = NULL,
                      maximum_number_of_iterations = 100L,
                      centroid_type = c("mean", "medoid"),
                      metric = c("l2", "pearson"),
                      warping_options = c(0.15, 0.15),
                      number_of_threads = 1L,
                      parallel_method = 0L,
                      distance_relative_tolerance = 0.001,
                      use_fence = FALSE,
                      check_total_dissimilarity = TRUE,
                      use_verbose = TRUE,
                      compute_overall_center = FALSE) {
  if (anyNA(x))
    cli::cli_abort("The input argument {.arg x} should not contain non-finite values.")

  if (anyNA(y))
    cli::cli_abort("The input argument {.arg y} should not contain non-finite values.")

  warping_class <- rlang::arg_match(warping_class)
  centroid_type <- rlang::arg_match(centroid_type)
  metric <- rlang::arg_match(metric)
  call <- rlang::call_match(defaults = TRUE)

  # Handle one-dimensional data
  if (length(dim(y)) == 2) {
    y <- array(y, c(dim(y)[1], 1, dim(y)[2]))
  }

  dims <- dim(y)
  N <- dims[1]
  L <- dims[2]
  M <- dims[3]

  # Handle vector grid
  if (is.vector(x)) {
    x <- matrix(x, N, M, byrow = TRUE)
  }

  # Handle seeds
  if (is.null(seeds)) {
    seeds <- sample(0:(N - 1), n_clusters)
  } else {
    seeds <- seeds - 1
  }

  # Compute common grid
  common_grid <- x[1, ]
  multiple_grids <- any(apply(x, 2, stats::sd) != 0)
  if (multiple_grids) {
    grid_min <- max(x[, 1])
    grid_max <- min(x[, M])
    common_grid <- seq(grid_min, grid_max, length.out = M)
  }

  if (warping_class == "srsf") {
    yperm <- aperm(y, c(2, 3, 1))

    if (multiple_grids) {
      for (l in 1:L) {
        for (n in 1:N)
          yperm[l, , n] <- stats::approx(
            x = x[n, ],
            y = yperm[l, , n],
            xout = common_grid
          )$y
      }
    }

    res <- fdasrvf::kmeans_align(
      f = yperm,
      time = common_grid,
      K = n_clusters,
      seeds = seeds + 1,
      centroid_type = centroid_type,
      max_iter = maximum_number_of_iterations,
      use_verbose = use_verbose
    )

    cluster_members <- lapply(1:n_clusters, function(k) which(res$labels == k))
    aligned_curves <- array(dim = c(N, L, M))
    warpings <- matrix(nrow = N, ncol = M)
    for (k in 1:n_clusters) {
      aligned_curves[cluster_members[[k]], , ] <- t(res$fn[[k]])
      warpings[cluster_members[[k]], ] <- t(res$gam[[k]])
    }

    out <- list(
      original_curves = aperm(array(res$f0, dim = c(L, M, N)), c(3, 1, 2)),
      aligned_curves = aligned_curves,
      center_curves = aperm(res$templates, c(3, 1, 2)),
      grid = res$time,
      n_clusters = n_clusters,
      memberships = res$labels,
      distances_to_center = res$distances_to_center,
      warpings = warpings,
      n_iterations = length(res$qun),
      call_name = rlang::call_name(call),
      call_args = rlang::call_args(call)
    )

    return(as_caps(out))
  }

  res <- kmap(
    x = x,
    y = y,
    seeds = seeds,
    warping_options = warping_options,
    n_clust = n_clusters,
    maximum_number_of_iterations = maximum_number_of_iterations,
    number_of_threads = number_of_threads,
    parallel_method = parallel_method,
    distance_relative_tolerance = distance_relative_tolerance,
    use_fence = use_fence,
    check_total_dissimilarity = check_total_dissimilarity,
    use_verbose = use_verbose,
    compute_overall_center = compute_overall_center,
    warping_method = warping_class,
    center_method = centroid_type,
    dissimilarity_method = metric,
    optimizer_method = "bobyqa"
  )

  original_curves <- res$y
  aligned_curves <- original_curves
  for (l in 1:L) {
    for (n in 1:N) {
      aligned_curves[n, l, ] <- stats::approx(
        x = res$x_final[n, ],
        y = original_curves[n, l, ],
        xout = common_grid
      )$y
      if (multiple_grids)
        original_curves[n, l, ] <- stats::approx(
          x = res$x[n, ],
          y = original_curves[n, l, ],
          xout = common_grid
        )$y
    }
  }

  centers <- res$y_centers_final
  for (l in 1:L) {
    for (k in 1:res$n_clust_final) {
      centers[k, l, ] <- stats::approx(
        x = res$x_centers_final[k, ],
        y = res$y_centers_final[k, l, ],
        xout = common_grid
      )$y
    }
  }

  warpings <- NULL
  if (warping_class == "none")
    warpings <- matrix(common_grid, nrow = N, ncol = M, byrow = TRUE)
  else {
    warpings <- matrix(nrow = N, ncol = M)
    for (n in 1:N) {
      if (warping_class == "shift")
        warpings[n, ] <- common_grid + res$parameters[n, 1]
      else if (warping_class == "dilation")
        warpings[n, ] <- common_grid * res$parameters[n, 1]
      else
        warpings[n, ] <- common_grid * res$parameters[n, 1] +
          res$parameters[n, 2]
    }
  }

  out <- list(
    original_curves = original_curves,
    aligned_curves = aligned_curves,
    center_curves = centers,
    grid = common_grid,
    n_clusters = res$n_clust_final,
    memberships = res$labels,
    distances_to_center = res$final_dissimilarity,
    warpings = warpings,
    n_iterations = res$iterations,
    call_name = rlang::call_name(call),
    call_args = rlang::call_args(call)
  )

  as_caps(out)
}

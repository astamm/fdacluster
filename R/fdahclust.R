#' Hierarchical Clustering for Functional Data with Amplitude and Phase Separation
#'
#' This function extends hierarchical agglomerative clustering to functional
#' data. It includes the possibility to separate amplitude and phase
#' information.
#'
#' The number of clusters is required as input because, with functional data,
#' once hierarchical clustering is performed, curves within clusters need to be
#' aligned to their corresponding centroid.
#'
#' @inheritParams fdakmeans
#' @param linkage_criterion A string specifying which linkage criterion should
#'   be used to compute distances between sets of curves. Choices are
#'   `"complete"` for complete linkage, `"average"` for average linkage and
#'   `"single"` for single linkage. See [stats::hclust()] for more details.
#'   Defqults to `"complete"`.
#'
#' @return An object of class [`caps`].
#'
#' @export
#' @examples
#' out <- fdahclust(simulated30$x, simulated30$y)
fdahclust <- function(x, y,
                      n_clusters = 1L,
                      warping_class = c("affine", "dilation", "none", "shift", "srsf"),
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
                      compute_overall_center = FALSE,
                      linkage_criterion = c("complete", "average", "single")) {
  if (anyNA(x))
    cli::cli_abort("The input argument {.arg x} should not contain non-finite values.")

  if (anyNA(y))
    cli::cli_abort("The input argument {.arg y} should not contain non-finite values.")

  warping_class <- rlang::arg_match(warping_class)
  centroid_type <- rlang::arg_match(centroid_type)
  metric <- rlang::arg_match(metric)
  linkage_criterion <- rlang::arg_match(linkage_criterion)
  call <- rlang::call_match(defaults = TRUE)

  # Handle one-dimensional data
  if (length(dim(y)) == 2) {
    y <- array(y, c(dim(y)[1], 1, dim(y)[2]))
  }

  dims <- dim(y)
  N <- dims[1]
  L <- dims[2]
  M <- dims[3]

  if (use_verbose)
    cli::cli_alert_info("Computing the distance matrix...")

  D <- fdadist(x = x, y = y, warping_class = warping_class, metric = metric)
  Dm <- as.matrix(D)

  if (use_verbose)
    cli::cli_alert_info("Calculating the tree...")

  hc <- stats::hclust(D, method = linkage_criterion)
  labels <- stats::cutree(tree = hc, k = n_clusters)

  if (use_verbose)
    cli::cli_alert_info("Aligning all curves with respect to their centroid...")

  kmresults <- lapply(1:n_clusters, function(k) {
    cluster_ids <- which(labels == k)
    medoid_idx <- which.min(rowSums(Dm[cluster_ids, cluster_ids]))
    fdakmeans(
      x = x[cluster_ids, ],
      y = y[cluster_ids, , ],
      n_clusters = 1,
      warping_class = warping_class,
      seeds = medoid_idx,
      maximum_number_of_iterations = maximum_number_of_iterations,
      centroid_type = centroid_type,
      metric = metric,
      warping_options = warping_options,
      number_of_threads = number_of_threads,
      parallel_method = parallel_method,
      distance_relative_tolerance = distance_relative_tolerance,
      use_fence = use_fence,
      check_total_dissimilarity = check_total_dissimilarity,
      use_verbose = FALSE,
      compute_overall_center = compute_overall_center
    )
  })

  if (use_verbose)
    cli::cli_alert_info("Consolidating output...")

  grid <- kmresults[[1]]$grid
  original_curves <- array(dim = c(N, L, M))
  aligned_curves <- array(dim = c(N, L, M))
  center_curves <- array(dim = c(n_clusters, L, M))
  warpings <- matrix(nrow = N, ncol = M)
  dtc <- numeric(N)
  for (k in 1:n_clusters) {
    cluster_ids <- which(labels == k)
    original_curves[cluster_ids, , ] <- kmresults[[k]]$original_curves
    aligned_curves[cluster_ids, , ] <- kmresults[[k]]$aligned_curves
    center_curves[k, , ] <- kmresults[[k]]$center_curves
    warpings[cluster_ids, ] <- kmresults[[k]]$warpings
    dtc[cluster_ids] <- kmresults[[k]]$distances_to_center
  }

  out <- list(
    original_curves = original_curves,
    aligned_curves = aligned_curves,
    center_curves = center_curves,
    grid = grid,
    n_clusters = n_clusters,
    memberships = labels,
    distances_to_center = dtc,
    warpings = warpings,
    n_iterations = 0,
    call_name = rlang::call_name(call),
    call_args = rlang::call_args(call)
  )

  return(as_caps(out))
}

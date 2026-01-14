#' Performs density-based clustering for functional data with amplitude and
#' phase separation
#'
#' This function extends `DBSCAN` to functional data. It includes the
#' possibility to separate amplitude and phase information.
#'
#' @inheritParams fdahclust
#'
#' @return An object of class [`caps`].
#'
#' @export
#' @examplesIf requireNamespace("lpSolve", quietly = TRUE)
#'
#' #----------------------------------
#' # Extracts 15 out of the 30 simulated curves in `simulated30_sub` data set
#' idx <- c(1:5, 11:15)
#' x <- simulated30_sub$x[idx, ]
#' y <- simulated30_sub$y[idx, , ]
#'
#' #----------------------------------
#' # Runs an HAC with affine alignment, searching for 2 clusters
#' out <- fdadbscan(
#'   x = x,
#'   y = y,
#'   warping_class = "affine",
#'   metric = "normalized_l2"
#' )
#'
#' #----------------------------------
#' # Then visualize the results
#' # Either with ggplot2 via ggplot2::autoplot(out)
#' # or using graphics::plot()
#' # You can visualize the original and aligned curves with:
#' plot(out, type = "amplitude")
#' # Or the estimated warping functions with:
#' plot(out, type = "phase")
fdadbscan <- function(
  x,
  y,
  is_domain_interval = FALSE,
  transformation = c("identity", "srvf"),
  warping_class = c("none", "shift", "dilation", "affine", "bpd"),
  centroid_type = "mean",
  metric = c("l2", "normalized_l2", "pearson"),
  cluster_on_phase = FALSE,
  use_verbose = FALSE,
  warping_options = c(0.15, 0.15),
  maximum_number_of_iterations = 100L,
  number_of_threads = 1L,
  parallel_method = 0L,
  distance_relative_tolerance = 0.001,
  use_fence = FALSE,
  check_total_dissimilarity = TRUE,
  compute_overall_center = FALSE
) {
  call <- rlang::call_match(defaults = TRUE)
  callname <- rlang::call_name(call)
  callargs <- rlang::call_args(call)

  transformation <- rlang::arg_match(transformation)
  callargs$transformation <- transformation

  warping_class <- rlang::arg_match(warping_class)
  callargs$warping_class <- warping_class

  metric <- rlang::arg_match(metric)
  callargs$metric <- metric

  l <- format_inputs(x, y, is_domain_interval)
  check_option_compatibility(
    is_domain_interval = is_domain_interval,
    transformation = transformation,
    warping_class = warping_class,
    metric = metric
  )

  x <- l$x
  y <- l$y
  dims <- dim(y)
  N <- dims[1]
  L <- dims[2]
  M <- dims[3]

  centroid_type_args <- check_centroid_type(centroid_type)
  centroid_name <- centroid_type_args$name
  centroid_extra <- centroid_type_args$extra

  if (centroid_name != "medoid" && parallel_method == 1L) {
    cli::cli_abort(
      "Parallelization on the distance calculation loop is only available for computing medoids."
    )
  }

  callargs$centroid_type <- centroid_name
  callargs$centroid_extra <- centroid_extra

  if (warping_class == "none" && cluster_on_phase) {
    cli::cli_abort(
      "It makes no sense to cluster based on phase variability if no alignment is performed."
    )
  }

  if (use_verbose) {
    cli::cli_alert_info("Computing the distance matrix...")
  }

  D <- fdadist(
    x = x,
    y = y,
    is_domain_interval = is_domain_interval,
    transformation = transformation,
    warping_class = warping_class,
    metric = metric,
    cluster_on_phase = cluster_on_phase
  )
  Dm <- as.matrix(D)

  results <- lapply(2:N, \(.min_pts) {
    dsts <- sort(dbscan::kNNdist(D, k = .min_pts - 1))
    eps <- dsts[which.max(diff(dsts))]
    obj <- dbscan::frNN(D, eps = eps)
    dbscan::dbscan(obj, minPts = .min_pts)
  })
  sils <- sapply(results, \(.res) {
    if (length(unique(.res$cluster)) == 1) {
      return(NA)
    }
    mean(cluster::silhouette(.res$cluster, D)[, "sil_width"])
  })

  if (all(is.na(sils))) {
    dbres <- results[[1]]
  } else {
    dbres <- results[[which.max(sils)]]
  }
  labels <- dbres$cluster
  n_clusters <- length(unique(labels[labels > 0]))

  if (use_verbose) {
    cli::cli_alert_info("Aligning all curves with respect to their centroid...")
  }

  kmresults <- lapply(1:n_clusters, function(k) {
    cluster_ids <- which(labels == k)
    medoid_idx <- which.min(rowSums(Dm[cluster_ids, cluster_ids, drop = FALSE]))
    fdakmeans(
      x = x[cluster_ids, , drop = FALSE],
      y = y[cluster_ids, , , drop = FALSE],
      n_clusters = 1,
      seeds = medoid_idx,
      is_domain_interval = is_domain_interval,
      transformation = transformation,
      warping_class = warping_class,
      centroid_type = centroid_type,
      metric = metric,
      maximum_number_of_iterations = maximum_number_of_iterations,
      warping_options = warping_options,
      number_of_threads = number_of_threads,
      parallel_method = parallel_method,
      distance_relative_tolerance = distance_relative_tolerance,
      cluster_on_phase = cluster_on_phase,
      use_fence = use_fence,
      check_total_dissimilarity = check_total_dissimilarity,
      use_verbose = FALSE,
      compute_overall_center = compute_overall_center,
      add_silhouettes = FALSE
    )
  })

  if (use_verbose) {
    cli::cli_alert_info("Consolidating output...")
  }

  original_curves <- array(dim = c(N, L, M))
  original_curves[labels == 0, , ] <- y[labels == 0, , ]
  original_grids <- matrix(nrow = N, ncol = M)
  original_grids[labels == 0, ] <- x[labels == 0, ]
  aligned_grids <- matrix(nrow = N, ncol = M)
  aligned_grids[labels == 0, ] <- x[labels == 0, ]
  center_curves <- array(dim = c(n_clusters, L, M))
  center_grids <- matrix(nrow = n_clusters, ncol = M)
  dtc <- numeric(N)
  dtc[labels == 0] <- 0
  for (k in 1:n_clusters) {
    cluster_ids <- which(labels == k)
    original_curves[cluster_ids, , ] <- kmresults[[k]]$original_curves
    original_grids[cluster_ids, ] <- kmresults[[k]]$original_grids
    aligned_grids[cluster_ids, ] <- kmresults[[k]]$aligned_grids
    center_curves[k, , ] <- kmresults[[k]]$center_curves
    center_grids[k, ] <- kmresults[[k]]$center_grids[1, ]
    dtc[cluster_ids] <- kmresults[[k]]$distances_to_center
  }

  silhouettes <- NULL
  if (n_clusters > 1) {
    D <- fdadist(
      x = aligned_grids,
      y = original_curves,
      is_domain_interval = is_domain_interval,
      transformation = transformation,
      warping_class = "none",
      metric = metric,
      cluster_on_phase = FALSE
    )
    silhouettes <- cluster::silhouette(labels, D)[, "sil_width"]
  }

  out <- list(
    original_curves = original_curves,
    original_grids = original_grids,
    aligned_grids = aligned_grids,
    center_curves = center_curves,
    center_grids = center_grids,
    n_clusters = n_clusters,
    memberships = labels,
    distances_to_center = dtc,
    silhouettes = silhouettes,
    amplitude_variation = sum(sapply(kmresults, \(.x) .x$amplitude_variation)),
    total_variation = sum(sapply(kmresults, \(.x) .x$total_variation)),
    n_iterations = 0,
    call_name = callname,
    call_args = callargs
  )

  as_caps(out)
}

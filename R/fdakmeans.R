#' Performs k-means clustering for functional data with amplitude and phase
#' separation
#'
#' This function provides implementations of the k-means clustering algorithm
#' for functional data, with possible joint amplitude and phase separation. A
#' number of warping class are implemented to achieve this separation.
#'
#' @param x A numeric vector of length \eqn{M} or a numeric matrix of shape
#'   \eqn{N \times M} or an object of class [`funData::funData`]. If a numeric
#'   vector or matrix, it specifies the grid(s) of size \eqn{M} on which each of
#'   the \eqn{N} curves have been observed. If an object of class
#'   [`funData::funData`], it contains the whole functional data set and the `y`
#'   argument is not used.
#' @param y Either a numeric matrix of shape \eqn{N \times M} or a numeric array
#'   of shape \eqn{N \times L \times M} or an object of class [`fda::fd`]. If a
#'   numeric matrix or array, it specifies the \eqn{N}-sample of
#'   \eqn{L}-dimensional curves observed on grids of size \eqn{M}. If an object
#'   of class [`fda::fd`], it contains all the necessary information about the
#'   functional data set to be able to evaluate it on user-defined grids.
#' @param n_clusters An integer value specifying the number of clusters.
#'   Defaults to `1L`.
#' @param seeds An integer value or vector specifying the indices of the initial
#'   centroids. If an integer vector, it is interpreted as the indices of the
#'   intial centroids and should therefore be of length `n_clusters`. If an
#'   integer value, it is interpreted as the index of the first initial centroid
#'   and subsequent centroids are chosen according to the k-means++ strategy. It
#'   can be `NULL` in which case the argument `seeding_strategy` is used to
#'   automatically provide suitable indices. Defaults to `NULL`.
#' @param seeding_strategy A character string specifying the strategy for
#'   choosing the initial centroids in case the argument `seeds` is set to
#'   `NULL`. Choices are
#'   [`"kmeans++"`](https://en.wikipedia.org/wiki/K-means%2B%2B),
#'   `"exhaustive-kmeans++"` which performs an exhaustive search over the choice
#'   of the first centroid, `"exhaustive"` which tries on all combinations of
#'   initial centroids or `"hclust"` which first performs hierarchical
#'   clustering using Ward's linkage criterion to identify initial centroids.
#'   Defaults to `"kmeans++"`, which is the fastest strategy.
#' @param warping_class A string specifying the warping class Choices are
#'   `"affine"`, `"dilation"`, `"none"`, `"shift"` or `"srsf"`. Defaults to
#'   `"affine"`. The SRSF class is the only class which is boundary-preserving.
#' @param centroid_type A string specifying the type of centroid to compute.
#'   Choices are `"mean"`, `"medoid"`, `"lowess"` or `"poly"`. Defaults to
#'   `"mean"`. If LOWESS appproximation is chosen, the user can append an
#'   integer between 0 and 100 as in `"lowess20"`. This number will be used as
#'   the smoother span. This gives the proportion of points in the plot which
#'   influence the smooth at each value. Larger values give more smoothness. The
#'   default value is 10%. If polynomial approximation is chosen, the user can
#'   append an positive integer as in `"poly3"`. This number will be used as the
#'   degree of the polynomial model. The default value is `4L`.
#' @param metric A string specifying the metric used to compare curves. Choices
#'   are `"l2"` or `"pearson"`. Defaults to `"l2"`. Used only when
#'   `warping_class != "srsf"`. For the boundary-preserving warping class, the
#'   L2 distance between the SRSFs of the original curves is used.
#' @param cluster_on_phase A boolean specifying whether clustering should be
#'   based on phase variation or amplitude variation. Defaults to `FALSE` which
#'   implies amplitude variation.
#' @param use_verbose A boolean specifying whether the algorithm should output
#'   details of the steps to the console. Defaults to `TRUE`.
#' @param warping_options A numeric vector supplied as a helper to the chosen
#'   `warping_class` to decide on warping parameter bounds. This is used only
#'   when `warping_class != "srsf"`.
#' @param maximum_number_of_iterations An integer specifying the maximum number
#'   of iterations before the algorithm stops if no other convergence criterion
#'   was met. Defaults to `100L`.
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
#' @param add_silhouettes A boolean specifying whether silhouette values should
#'   be computed for each observation for internal validation of the clustering
#'   structure. Defaults to `TRUE`.
#' @param expand_domain A boolean specifying how to define the within-cluster
#'   common grids. When set to `FALSE`, the intersection of the individual
#'   domains is used. When set to `TRUE`, the union is used and mean imputation
#'   is performed to fill in the missing values. Defaults to `TRUE`.
#'
#' @return An object of class [`caps`].
#'
#' @export
#' @examples
#' #----------------------------------
#' # Extracts 15 out of the 30 simulated curves in `simulated30_sub` data set
#' idx <- c(1:5, 11:15, 21:25)
#' x <- simulated30_sub$x[idx, ]
#' y <- simulated30_sub$y[idx, , ]
#'
#' #----------------------------------
#' # Runs a k-means clustering with affine alignment, searching for 2 clusters
#' out <- fdakmeans(
#'   x = x,
#'   y = y,
#'   n_clusters = 2,
#'   warping_class = "affine"
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
fdakmeans <- function(x, y = NULL,
                      n_clusters = 1L,
                      seeds = NULL,
                      seeding_strategy = c("kmeans++", "exhaustive-kmeans++", "exhaustive", "hclust"),
                      warping_class = c("affine", "dilation", "none", "shift", "srsf"),
                      centroid_type = "mean",
                      metric = c("l2", "pearson"),
                      cluster_on_phase = FALSE,
                      use_verbose = TRUE,
                      warping_options = c(0.15, 0.15),
                      maximum_number_of_iterations = 100L,
                      number_of_threads = 1L,
                      parallel_method = 0L,
                      distance_relative_tolerance = 0.001,
                      use_fence = FALSE,
                      check_total_dissimilarity = TRUE,
                      compute_overall_center = FALSE,
                      add_silhouettes = TRUE,
                      expand_domain = TRUE) {
  call <- rlang::call_match(defaults = TRUE)
  call_name <- rlang::call_name(call)
  call_args <- rlang::call_args(call)

  l <- format_inputs(x, y)
  x <- l$x
  y <- l$y
  dims <- dim(y)
  N <- dims[1]
  L <- dims[2]
  M <- dims[3]

  seeding_strategy <- rlang::arg_match(seeding_strategy)
  warping_class <- rlang::arg_match(warping_class)
  metric <- rlang::arg_match(metric)

  centroid_type_args <- check_centroid_type(centroid_type)
  centroid_name <- centroid_type_args$name
  centroid_extra <- centroid_type_args$extra

  if (centroid_name != "medoid" && parallel_method == 1L)
    cli::cli_abort("Parallelization on the distance calculation loop is only available for computing medoids.")

  if (warping_class == "none" && cluster_on_phase)
    cli::cli_abort("It makes no sense to cluster based on phase variability if no alignment is performed.")

  # Handle seeds
  if (is.null(seeds)) {
    if (use_verbose)
      cli::cli_alert_info("Computing initial centroids using {seeding_strategy} strategy...")
    if (seeding_strategy == "hclust") {
      out <- fdahclust(
        x = x,
        y = y,
        n_clusters = n_clusters,
        warping_class = warping_class,
        maximum_number_of_iterations = maximum_number_of_iterations,
        centroid_type = centroid_type,
        metric = metric,
        linkage_criterion = "ward.D2",
        cluster_on_phase = cluster_on_phase,
        use_verbose = FALSE
      )
      seeds <- 1:n_clusters |>
        purrr::map(~ which(out$memberships == .x)) |>
        purrr::map_int(~ .x[which.min(out$distances_to_center[.x])])
    } else if (seeding_strategy == "kmeans++") {
      D <- fdadist(
        x = x,
        y = y,
        warping_class = warping_class,
        metric = metric,
        cluster_on_phase = cluster_on_phase
      )
      Dm <- as.matrix(D)
      seeds <- sample(1:N, 1L)
      if (n_clusters > 1L) {
        for (k in 2:n_clusters) {
          Dsub <- Dm[seeds, -seeds, drop = FALSE]
          Dvec <- apply(Dsub, 2L, min)
          non_seeds <- setdiff(1:N, seeds)
          seeds <- c(seeds, sample(non_seeds, 1L, prob = Dvec^2))
        }
      }
    } else if (seeding_strategy == "exhaustive-kmeans++") {
      D <- fdadist(
        x = x,
        y = y,
        warping_class = warping_class,
        metric = metric,
        cluster_on_phase = cluster_on_phase
      )
      Dm <- as.matrix(D)
      pb <- progressr::progressor(steps = N)
      out <- furrr::future_map(1:N, \(n) {
        pb()
        seeds <- n
        if (n_clusters > 1L) {
          for (k in 2:n_clusters) {
            Dsub <- Dm[seeds, -seeds, drop = FALSE]
            Dvec <- apply(Dsub, 2L, min)
            non_seeds <- setdiff(1:N, seeds)
            seeds <- c(seeds, sample(non_seeds, 1L, prob = Dvec^2))
          }
        }
        km <- fdakmeans(
          x = x,
          y = y,
          n_clusters = n_clusters,
          warping_class = warping_class,
          seeds = seeds,
          maximum_number_of_iterations = maximum_number_of_iterations,
          centroid_type = centroid_type,
          metric = metric,
          warping_options = warping_options,
          distance_relative_tolerance = distance_relative_tolerance,
          use_fence = use_fence,
          cluster_on_phase = cluster_on_phase,
          use_verbose = FALSE,
          add_silhouettes = add_silhouettes
        )
        list(caps = km, totss = sum(km$distances_to_center))
      }, .options = furrr::furrr_options(seed = TRUE))
      best_idx <- which.min(purrr::map_dbl(out, "totss"))
      return(purrr::map(out, "caps")[[best_idx]])
    } else if (seeding_strategy == "exhaustive") {
      sols <- utils::combn(N, n_clusters, simplify = FALSE)
      pb <- progressr::progressor(steps = length(sols))
      sols <- furrr::future_map(sols, \(.seeds) {
        pb()
        fdakmeans(
          x = x, y = y,
          n_clusters = n_clusters,
          warping_class = warping_class,
          seeds = .seeds,
          maximum_number_of_iterations = maximum_number_of_iterations,
          centroid_type = centroid_type,
          metric = metric,
          warping_options = warping_options,
          distance_relative_tolerance = distance_relative_tolerance,
          cluster_on_phase = cluster_on_phase,
          use_fence = use_fence,
          check_total_dissimilarity = check_total_dissimilarity,
          use_verbose = FALSE,
          add_silhouettes = FALSE
        )
      }, .options = furrr::furrr_options(seed = NULL))
      dtcs <- sols |>
        purrr::map("distances_to_center") |>
        purrr::map_dbl(sum)
      return(as_caps(sols[[which.min(dtcs)]]))
    }
  } else {
    n_centroids <- length(seeds)
    if (n_centroids != n_clusters && n_centroids != 1L)
      cli::cli_abort("The number of initial centroid indices provided by the {.arg seeds} argument should be either 1 or {n_clusters}.")
    if (n_centroids == 1L && n_clusters > 1L) {
      D <- fdadist(
        x = x,
        y = y,
        warping_class = warping_class,
        metric = metric,
        cluster_on_phase = cluster_on_phase
      )
      Dm <- as.matrix(D)
      for (k in 2:n_clusters) {
        Dsub <- Dm[seeds, -seeds, drop = FALSE]
        Dvec <- apply(Dsub, 2L, min)
        non_seeds <- setdiff(1:N, seeds)
        seeds <- c(seeds, sample(non_seeds, 1L, prob = Dvec^2))
      }
    }
  }

  call_args$seeds <- seeds
  seeds <- seeds - 1

  if (warping_class == "srsf") {
    # Compute common grid
    common_grid <- x[1, ]
    multiple_grids <- FALSE
    if (N > 1) {
      multiple_grids <- any(apply(x, 2, stats::sd) != 0)
      if (multiple_grids) {
        grid_min <- max(x[, 1])
        grid_max <- min(x[, M])
        common_grid <- seq(grid_min, grid_max, length.out = M)
      }
    }

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

    original_curves <- aperm(array(res$f0, dim = c(L, M, N)), c(3, 1, 2))
    cluster_members <- lapply(1:n_clusters, function(k) which(res$labels == k))
    aligned_curves <- array(dim = c(N, L, M))
    center_curves <- aperm(res$templates, c(3, 1, 2))
    warpings <- matrix(nrow = N, ncol = M)
    for (k in 1:n_clusters) {
      aligned_curves[cluster_members[[k]], , ] <- t(res$fn[[k]])
      warpings[cluster_members[[k]], ] <- t(res$gam[[k]])
    }

    q0 <- res$q0
    if (length(dim(q0)) == 2L) # This should be done in fdasrvf package
      dim(q0) <- c(1, dim(q0))
    amplitude_variation <- sum(res$distances_to_center^2)
    total_variation <- sum(purrr::map_dbl(1:N, \(n) {
      sum(purrr::map_dbl(1:L, \(l) {
        trapz(
          x = common_grid,
          y = (q0[l, , n] - res$templates.q[l, , res$labels[n]])^2
        )
      }))
    }))

    silhouettes <- NULL
    if (n_clusters > 1 && add_silhouettes) {
      D <- fdadist(
        x = common_grid,
        y = aligned_curves,
        warping_class = "none",
        metric = metric
      )
      silhouettes <- cluster::silhouette(res$labels, D)[, "sil_width"]
    }

    out <- list(
      original_curves = original_curves,
      aligned_curves = aligned_curves,
      center_curves = center_curves,
      warpings = warpings,
      grids = matrix(res$time, nrow = n_clusters, ncol = M, byrow = TRUE),
      n_clusters = n_clusters,
      memberships = res$labels,
      distances_to_center = res$distances_to_center,
      silhouettes = silhouettes,
      amplitude_variation = amplitude_variation,
      total_variation = total_variation,
      n_iterations = length(res$qun),
      call_name = call_name,
      call_args = call_args
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
    center_args = centroid_extra,
    cluster_on_phase = cluster_on_phase,
    use_fence = use_fence,
    check_total_dissimilarity = check_total_dissimilarity,
    use_verbose = use_verbose,
    compute_overall_center = compute_overall_center,
    warping_method = warping_class,
    center_method = centroid_name,
    dissimilarity_method = metric,
    optimizer_method = "bobyqa"
  )

  # Compute common grid per cluster
  common_grids <- purrr::map(1:n_clusters, \(cluster_id) {
    grids <- res$x_final[res$labels == cluster_id, , drop = FALSE]
    common_grid <- grids[1, ]
    if (nrow(grids) > 1) {
      multiple_grids <- any(apply(grids, 2, stats::sd) != 0)
      if (multiple_grids) {
        if (expand_domain) {
          grid_min <- min(grids[, 1])
          grid_max <- max(grids[, M])
        } else {
          grid_min <- max(grids[, 1])
          grid_max <- min(grids[, M])
        }
        common_grid <- seq(grid_min, grid_max, length.out = M)
      }
    }
    common_grid
  })
  common_grids <- do.call(rbind, common_grids)

  original_curves <- res$y
  aligned_curves <- original_curves
  for (l in 1:L) {
    for (n in 1:N) {
      aligned_curves[n, l, ] <- stats::approx(
        x = res$x_final[n, ],
        y = original_curves[n, l, ],
        xout = common_grids[res$labels[n], ]
      )$y
      original_curves[n, l, ] <- stats::approx(
        x = res$x[n, ],
        y = original_curves[n, l, ],
        xout = common_grids[res$labels[n], ]
      )$y
    }
  }

  centers <- res$y_centers_final
  for (l in 1:L) {
    for (k in 1:res$n_clust_final) {
      centers[k, l, ] <- stats::approx(
        x = res$x_centers_final[k, ],
        y = res$y_centers_final[k, l, ],
        xout = common_grids[k, ]
      )$y
    }
  }

  if (expand_domain) {
    # Mean imputation
    for (l in 1:L) {
      X <- aligned_curves[, l, ]
      Xc <- centers[, l, ]
      if (n_clusters == 1) Xc <- matrix(Xc, nrow = 1)
      X <- rbind(X, Xc)
      Xi <- impute_via_mean(X, c(res$labels, 1:n_clusters))
      aligned_curves[, l, ] <- Xi[1:N, ]
      centers[, l, ] <- Xi[(N + 1):nrow(Xi), ]
    }
  }

  warpings <- matrix(nrow = N, ncol = M)
  if (warping_class == "none") {
    for (n in 1:N)
      warpings[n, ] <- common_grids[res$labels[n], ]
  } else {
    for (n in 1:N) {
      if (warping_class == "shift")
        warpings[n, ] <- common_grids[res$labels[n], ] + res$parameters[n, 1]
      else if (warping_class == "dilation")
        warpings[n, ] <- common_grids[res$labels[n], ] * res$parameters[n, 1]
      else
        warpings[n, ] <- common_grids[res$labels[n], ] * res$parameters[n, 1] +
          res$parameters[n, 2]
    }
  }

  silhouettes <- NULL
  if (n_clusters > 1 && add_silhouettes) {
    D <- fdadist(
      x = res$x_final,
      y = res$y,
      warping_class = "none",
      metric = metric
    )
    silhouettes <- cluster::silhouette(res$labels, D)[, "sil_width"]
  }

  out <- list(
    original_curves = original_curves,
    aligned_curves = aligned_curves,
    center_curves = centers,
    warpings = warpings,
    grids = common_grids,
    n_clusters = res$n_clust_final,
    memberships = res$labels,
    distances_to_center = res$final_dissimilarity,
    silhouettes = silhouettes,
    amplitude_variation = res$amplitude_variation,
    total_variation = res$total_variation,
    n_iterations = res$iterations,
    call_name = call_name,
    call_args = call_args
  )

  as_caps(out)
}

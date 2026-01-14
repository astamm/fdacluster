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
#' @param is_domain_interval A boolean specifying whether the sample of curves
#'   is defined on a fixed interval. Defaults to `FALSE`.
#' @param transformation A string specifying the transformation to apply to the
#'   original sample of curves. Choices are no transformation (`transformation =
#'   "identity"`) or square-root velocity function `transformation = "srvf"`.
#'   Defaults to `"identity"`.
#' @param warping_class A string specifying the class of warping functions.
#'   Choices are no warping (`warping_class = "none"`), shift `y = x + b`
#'   (`warping_class = "shift"`), dilation `y = ax` (`warping_class =
#'   "dilation"`), affine `y = ax + b` (`warping_class = "affine"`) or
#'   boundary-preserving diffeomorphism (`warping_class = "bpd"`). Defaults to
#'   `"none"`.
#' @param centroid_type A string specifying the type of centroid to compute.
#'   Choices are `"mean"`, `"median"` `"medoid"`, `"lowess"` or `"poly"`.
#'   Defaults to `"mean"`. If LOWESS appproximation is chosen, the user can
#'   append an integer between 0 and 100 as in `"lowess20"`. This number will be
#'   used as the smoother span. This gives the proportion of points in the plot
#'   which influence the smooth at each value. Larger values give more
#'   smoothness. The default value is 10%. If polynomial approximation is
#'   chosen, the user can append an positive integer as in `"poly3"`. This
#'   number will be used as the degree of the polynomial model. The default
#'   value is `4L`.
#' @param metric A string specifying the metric used to compare curves. Choices
#'   are `"l2"`, `"normalized_l2"` or `"pearson"`. If `transformation ==
#'   "srvf"`, the metric **must be** `"l2"` because the SRVF transform maps
#'   absolutely continuous functions to square-integrable functions. If
#'   `transformation == "identity"` and `warping_class` is either `dilation` or
#'   `affine`, the metric cab be either `"normalized_l2"` or `"pearson"`. The L2
#'   distance is indeed **not** dilation-invariant or affine-invariant. The
#'   metric can also be `"l2"` if `warping_class == "shift"`. Defaults to
#'   `"l2"`.
#' @param cluster_on_phase A boolean specifying whether clustering should be
#'   based on phase variation or amplitude variation. Defaults to `FALSE` which
#'   implies amplitude variation.
#' @param use_verbose A boolean specifying whether the algorithm should output
#'   details of the steps to the console. Defaults to `FALSE`.
#' @param warping_options A numeric vector supplied as a helper to the chosen
#'   `warping_class` to decide on warping parameter bounds. This is used only
#'   when `warping_class != "srvf"`.
#' @param maximum_number_of_iterations An integer specifying the maximum number
#'   of iterations before the algorithm stops if no other convergence criterion
#'   was met. Defaults to `100L`.
#' @param number_of_threads An integer value specifying the number of threads
#'   used for parallelization. Defaults to `1L`. This is used only when
#'   `warping_class != "srvf"`.
#' @param parallel_method An integer value specifying the type of desired
#'   parallelization for template computation, If `0L`, templates are computed
#'   in parallel. If `1L`, parallelization occurs within a single template
#'   computation (only for the medoid method as of now). Defaults to `0L`. This
#'   is used only when `warping_class != "srvf"`.
#' @param distance_relative_tolerance A numeric value specifying a relative
#'   tolerance on the distance update between two iterations. If all
#'   observations have not sufficiently improved in that sense, the algorithm
#'   stops. Defaults to `1e-3`. This is used only when `warping_class !=
#'   "srvf"`.
#' @param use_fence A boolean specifying whether the fence algorithm should be
#'   used to robustify the algorithm against outliers. Defaults to `FALSE`. This
#'   is used only when `warping_class != "srvf"`.
#' @param check_total_dissimilarity A boolean specifying whether an additional
#'   stopping criterion based on improvement of the total dissimilarity should
#'   be used. Defaults to `TRUE`. This is used only when `warping_class !=
#'   "srvf"`.
#' @param compute_overall_center A boolean specifying whether the overall center
#'   should be also computed. Defaults to `FALSE`. This is used only when
#'   `warping_class != "srvf"`.
#' @param add_silhouettes A boolean specifying whether silhouette values should
#'   be computed for each observation for internal validation of the clustering
#'   structure. Defaults to `TRUE`.
#'
#' @return An object of class [`caps`].
#'
#' @export
#' @examplesIf requireNamespace("lpSolve", quietly = TRUE)
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
fdakmeans <- function(
  x,
  y = NULL,
  n_clusters = 1L,
  seeds = NULL,
  seeding_strategy = c(
    "kmeans++",
    "exhaustive-kmeans++",
    "exhaustive",
    "hclust"
  ),
  is_domain_interval = FALSE,
  transformation = c("identity", "srvf"),
  warping_class = c("none", "shift", "dilation", "affine", "bpd"),
  centroid_type = "mean",
  metric = c("l2", "normalized_l2", "pearson"),
  cluster_on_phase = FALSE, # TODO: implement arg to say if distance is based on original data, amplitude or phase
  use_verbose = FALSE,
  warping_options = c(0.15, 0.15),
  maximum_number_of_iterations = 100L,
  number_of_threads = 1L,
  parallel_method = 0L,
  distance_relative_tolerance = 0.001,
  use_fence = FALSE,
  check_total_dissimilarity = TRUE,
  compute_overall_center = FALSE,
  add_silhouettes = TRUE
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

  seeding_strategy <- rlang::arg_match(seeding_strategy)
  callargs$seeding_strategy <- seeding_strategy

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

  # Handle seeds
  if (is.null(seeds)) {
    if (use_verbose) {
      cli::cli_alert_info(
        "Computing initial centroids using {seeding_strategy} strategy..."
      )
    }
    if (seeding_strategy == "hclust") {
      out <- fdahclust(
        x = x,
        y = y,
        n_clusters = n_clusters,
        is_domain_interval = is_domain_interval,
        transformation = transformation,
        warping_class = warping_class,
        metric = metric,
        cluster_on_phase = cluster_on_phase,
        maximum_number_of_iterations = maximum_number_of_iterations,
        centroid_type = centroid_type,
        linkage_criterion = "ward.D2",
        use_verbose = FALSE
      )
      seeds <- 1:n_clusters |>
        lapply(\(.x) which(out$memberships == .x)) |>
        sapply(\(.x) .x[which.min(out$distances_to_center[.x])])
    } else if (seeding_strategy == "kmeans++") {
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
        is_domain_interval = is_domain_interval,
        transformation = transformation,
        warping_class = warping_class,
        metric = metric,
        cluster_on_phase = cluster_on_phase
      )
      Dm <- as.matrix(D)
      pb <- progressr::progressor(steps = N)
      out <- future.apply::future_lapply(
        1:N,
        \(n) {
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
            is_domain_interval = is_domain_interval,
            transformation = transformation,
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
            add_silhouettes = FALSE
          )
          list(caps = km, totss = sum(km$distances_to_center))
        },
        future.seed = TRUE,
        future.packages = "fdacluster"
      )
      best_idx <- which.min(sapply(out, \(.x) .x$totss))
      return(lapply(out, \(.x) .x$caps)[[best_idx]])
    } else if (seeding_strategy == "exhaustive") {
      sols <- utils::combn(N, n_clusters, simplify = FALSE)
      pb <- progressr::progressor(steps = length(sols))
      sols <- future.apply::future_lapply(
        sols,
        \(.seeds) {
          pb()
          fdakmeans(
            x = x,
            y = y,
            n_clusters = n_clusters,
            is_domain_interval = is_domain_interval,
            transformation = transformation,
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
        },
        future.seed = TRUE,
        future.packages = "fdacluster"
      )
      dtcs <- sols |>
        lapply(\(.x) .x$distances_to_center) |>
        sapply(sum)
      return(as_caps(sols[[which.min(dtcs)]]))
    }
  } else {
    n_centroids <- length(seeds)
    if (n_centroids != n_clusters && n_centroids != 1L) {
      cli::cli_abort(
        "The number of initial centroid indices provided by the {.arg seeds} argument should be either 1 or {n_clusters}."
      )
    }
    if (n_centroids == 1L && n_clusters > 1L) {
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
      for (k in 2:n_clusters) {
        Dsub <- Dm[seeds, -seeds, drop = FALSE]
        Dvec <- apply(Dsub, 2L, min)
        non_seeds <- setdiff(1:N, seeds)
        seeds <- c(seeds, sample(non_seeds, 1L, prob = Dvec^2))
      }
    }
  }

  callargs$seeds <- seeds
  seeds <- seeds - 1

  if (transformation == "srvf" && warping_class %in% c("none", "bpd")) {
    if (!(centroid_name %in% c("mean", "medoid"))) {
      cli::cli_abort(
        "Only mean and medoid centroids are available for SRVFs using the {.fn fdasrvf::kmeans_align} function."
      )
    }

    yperm <- aperm(y, c(2, 3, 1))
    common_grid <- x[1, ]

    # AST: nonempty should be 1 but badly handled in fdasrvf
    res <- fdasrvf::kmeans_align(
      f = yperm,
      time = common_grid,
      K = n_clusters,
      seeds = seeds + 1,
      centroid_type = centroid_type,
      max_iter = maximum_number_of_iterations,
      use_verbose = use_verbose,
      nonempty = 2L,
      scale = FALSE,
      alignment = (warping_class == "bpd")
    )

    original_curves <- aperm(array(res$f0, dim = c(L, M, N)), c(3, 1, 2))
    cluster_members <- lapply(1:n_clusters, function(k) which(res$labels == k))
    aligned_curves <- array(dim = c(N, L, M))
    center_curves <- aperm(res$templates, c(3, 1, 2))
    warpings <- matrix(nrow = N, ncol = M)
    for (k in 1:n_clusters) {
      Nc <- length(cluster_members[[k]])
      aligned_curves[cluster_members[[k]], , ] <- aperm(
        array(res$fn[[k]], dim = c(L, M, Nc)),
        c(3, 1, 2)
      )
      warpings[cluster_members[[k]], ] <- t(apply(
        t(res$gam[[k]]),
        1,
        fdasrvf::invertGamma
      )) *
        diff(range(common_grid)) +
        min(common_grid)
    }

    if (cluster_on_phase) {
      res$distances_to_center <- sqrt(sapply(1:N, \(n) {
        trapz(
          x = common_grid,
          y = (warpings[n, ] - common_grid)^2
        )
      }))
    }

    # AST: warning -- if cluster on phase in ON, below does not make sense, also
    # clustering is still performed on amplitude for SRVF, only distances to
    # centers are transformed in the output
    q0 <- res$q0
    if (length(dim(q0)) == 2L) {
      # This should be done in fdasrvf package
      dim(q0) <- c(1, dim(q0))
    }
    amplitude_variation <- sum(res$distances_to_center^2)
    total_variation <- sum(sapply(1:N, \(n) {
      sum(sapply(1:L, \(l) {
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
        is_domain_interval = is_domain_interval,
        transformation = transformation,
        warping_class = "none",
        metric = metric,
        cluster_on_phase = FALSE
      )
      silhouettes <- cluster::silhouette(res$labels, D)[, "sil_width"]
    }

    out <- list(
      original_curves = original_curves,
      original_grids = matrix(res$time, nrow = N, ncol = M, byrow = TRUE),
      aligned_grids = warpings,
      center_curves = center_curves,
      center_grids = matrix(
        res$time,
        nrow = n_clusters,
        ncol = M,
        byrow = TRUE
      ),
      n_clusters = n_clusters,
      memberships = res$labels,
      distances_to_center = res$distances_to_center,
      silhouettes = silhouettes,
      amplitude_variation = amplitude_variation,
      total_variation = total_variation,
      n_iterations = length(res$qun),
      call_name = callname,
      call_args = callargs
    )

    return(as_caps(out))
  }

  # Handle transformation
  if (transformation == "srvf") {
    for (i in 1:N) {
      y[i, , ] <- fdasrvf::f_to_srvf(
        f = y[i, , ],
        time = x[i, ],
        multidimensional = TRUE
      )
    }
  }

  res <- kmap(
    x = x,
    y = y,
    n_clust = n_clusters,
    seeds = seeds,
    cluster_on_phase = cluster_on_phase,
    warping_options = warping_options,
    maximum_number_of_iterations = maximum_number_of_iterations,
    number_of_threads = number_of_threads,
    parallel_method = parallel_method,
    distance_relative_tolerance = distance_relative_tolerance,
    center_args = centroid_extra,
    use_fence = use_fence,
    check_total_dissimilarity = check_total_dissimilarity,
    use_verbose = use_verbose,
    compute_overall_center = compute_overall_center,
    warping_method = warping_class,
    center_method = centroid_name,
    dissimilarity_method = metric,
    optimizer_method = "bobyqa"
  )

  silhouettes <- NULL
  if (n_clusters > 1 && add_silhouettes) {
    D <- fdadist(
      x = res$x_final,
      y = res$y,
      is_domain_interval = is_domain_interval,
      transformation = transformation,
      warping_class = "none",
      metric = metric,
      cluster_on_phase = FALSE
    )
    silhouettes <- cluster::silhouette(res$labels, D)[, "sil_width"]
  }

  out <- list(
    original_curves = res$y,
    original_grids = res$x,
    aligned_grids = res$x_final,
    center_curves = res$y_centers_final,
    center_grids = res$x_centers_final,
    n_clusters = res$n_clust_final,
    memberships = res$labels,
    distances_to_center = res$final_dissimilarity,
    silhouettes = silhouettes,
    amplitude_variation = res$amplitude_variation,
    total_variation = res$total_variation,
    n_iterations = res$iterations,
    call_name = callname,
    call_args = callargs
  )

  as_caps(out)
}

#' Computes the distance matrix for functional data with amplitude and phase
#' separation
#'
#' This function computes the matrix of pairwise distances between curves a
#' functional data sample. This can be achieved with or without phase and
#' amplitude separation, which can be done using a variety of warping classes.
#'
#' @inheritParams fdakmeans
#' @param labels A character vector specifying curve labels. Defaults to `NULL`
#'   which uses sequential numbers as labels.
#'
#' @return A [stats::dist] object storing the distance matrix between the input
#'   curves using the metric specified through the argument `metric` and the
#'   warping class specified by the argument `warping_class`.
#'
#' @export
#' @examples
#' idx <- c(1:5, 11:15, 21:25)
#' D <- fdadist(simulated30_sub$x[idx, ], simulated30_sub$y[idx, , ])
fdadist <- function(x, y = NULL,
                    is_domain_interval = FALSE,
                    transformation = c("identity", "srvf"),
                    warping_class = c("none", "shift", "dilation", "affine", "bpd"),
                    metric = c("l2", "normalized_l2", "pearson"),
                    cluster_on_phase = FALSE,
                    labels = NULL) {
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

  if (warping_class == "none" && cluster_on_phase)
    cli::cli_abort("It makes no sense to cluster based on phase variability if no alignment is performed.")

  if (is.null(labels))
    labels <- 1:N

  index_table <- linear_index(N)
  curve_pair <- array(dim = c(2, L, M))
  grid_pair <- array(dim = c(2, M))

  .pairwise_distances <- function(index_table) {
    pb <- progressr::progressor(steps = nrow(index_table))
    future_map2_dbl(index_table$i, index_table$j, \(i, j) {
      pb()

      curve_pair[1, , ] <- y[i, , ]
      curve_pair[2, , ] <- y[j, , ]
      grid_pair[1, ] <- x[i, ]
      grid_pair[2, ] <- x[j, ]

      km <- fdakmeans(
        grid_pair,
        curve_pair,
        n_clusters = 1,
        seeds = 1,
        is_domain_interval = is_domain_interval,
        transformation = transformation,
        warping_class = warping_class,
        centroid_type = "medoid",
        metric = metric,
        cluster_on_phase = cluster_on_phase,
        use_verbose = FALSE,
        add_silhouettes = FALSE
      )

      out <- max(km$distances_to_center)

      km <- fdakmeans(
        grid_pair,
        curve_pair,
        seeds = 2,
        n_clusters = 1,
        is_domain_interval = is_domain_interval,
        transformation = transformation,
        warping_class = warping_class,
        centroid_type = "medoid",
        metric = metric,
        cluster_on_phase = cluster_on_phase,
        use_verbose = FALSE,
        add_silhouettes = FALSE
      )

      if (max(km$distances_to_center) < out)
        out <- max(km$distances_to_center)

      out
    }, future_seed = TRUE, future_packages = "fdacluster"
  )}

  d <- .pairwise_distances(index_table)

  attributes(d) <- NULL
  attr(d, "Labels") <- labels
  attr(d, "Size") <- N
  attr(d, "Diag") <- FALSE
  attr(d, "Upper") <- FALSE
  attr(d, "call") <- call
  attr(d, "method") <- metric
  class(d) <- "dist"
  d
}

linear_index <- function(n) {
  df <- expand.grid(j = 1:n, i = 1:n)
  df <- subset(df, df$j > df$i)
  rownames(df) <- NULL
  df[, 2:1]
}

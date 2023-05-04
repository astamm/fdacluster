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
                    warping_class = c("affine", "dilation", "none", "shift", "srsf"),
                    metric = c("l2", "pearson"),
                    cluster_on_phase = FALSE,
                    labels = NULL) {
  call <- rlang::call_match(defaults = TRUE)

  l <- format_inputs(x, y)
  x <- l$x
  y <- l$y
  dims <- dim(y)
  N <- dims[1]
  L <- dims[2]
  M <- dims[3]

  warping_class <- rlang::arg_match(warping_class)
  metric <- rlang::arg_match(metric)

  if (warping_class == "none" && cluster_on_phase)
    cli::cli_abort("It makes no sense to cluster based on phase variability if no alignment is performed.")

  if (is.null(labels))
    labels <- 1:N

  index_table <- linear_index(N)
  curve_pair <- array(dim = c(2, L, M))
  grid_pair <- array(dim = c(2, M))

  .pairwise_distances <- function(index_table) {
    pb <- progressr::progressor(steps = nrow(index_table))
    furrr::future_map2_dbl(index_table$i, index_table$j, \(i, j) {
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
        centroid_type = "medoid",
        warping_class = warping_class,
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
        centroid_type = "medoid",
        warping_class = warping_class,
        metric = metric,
        cluster_on_phase = cluster_on_phase,
        use_verbose = FALSE,
        add_silhouettes = FALSE
      )

      if (max(km$distances_to_center) < out)
        out <- max(km$distances_to_center)

      out
    }, .options = furrr::furrr_options(seed = TRUE))
  }

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
  res <- tidyr::expand_grid(i = 1:n, j = 1:n)
  res <- subset(res, res$j > res$i)
  # res$k <- as.integer(n * (res$i - 1) - res$i * (res$i - 1) / 2 + res$j - res$i)
  res
}

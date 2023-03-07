#' Distance Matrix for Functional Data with Amplitude and Phase Separation
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
#' D <- fdadist(simulated30$x, simulated30$y)
fdadist <- function(x, y,
                    warping_class = c("affine", "dilation", "none", "shift", "srsf"),
                    metric = c("l2", "pearson"),
                    labels = NULL) {
  if (anyNA(x))
    cli::cli_abort("The input argument {.arg x} should not contain non-finite values.")

  if (anyNA(y))
    cli::cli_abort("The input argument {.arg y} should not contain non-finite values.")

  warping_class <- rlang::arg_match(warping_class)
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

  if (is.null(labels))
    labels <- 1:N

  indices <- linear_index(N)
  curve_pair <- array(dim = c(2, L, M))
  grid_pair <- array(dim = c(2, M))

  .pairwise_distances <- function(linear_indices) {
    pb <- progressr::progressor(along = linear_indices)
    furrr::future_map_dbl(linear_indices, ~ {
      pb()

      i <- indices$i[.x]
      j <- indices$j[.x]
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
        use_verbose = FALSE,
        add_silhouettes = FALSE
      )

      if (max(km$distances_to_center) < out)
        out <- max(km$distances_to_center)

      out
    }, .options = furrr::furrr_options(seed = TRUE))
  }

  d <- .pairwise_distances(indices$k)

  attributes(d) <- NULL
  attr(d, "Labels") <- labels
  attr(d, "Size") <- N
  attr(d, "Diag") <- FALSE
  attr(d, "Upper") <- FALSE
  attr(d, "call") <- rlang::call_match(defaults = TRUE)
  attr(d, "method") <- metric
  class(d) <- "dist"
  d
}

linear_index <- function(n) {
  res <- tidyr::expand_grid(i = 1:n, j = 1:n)
  res <- subset(res, res$j > res$i)
  res$k <- n * (res$i - 1) - res$i * (res$i - 1) / 2 + res$j - res$i
  res
}

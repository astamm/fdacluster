#' Multiple CAPS `mcaps` objects
#'
#' This function searches for clusters in the input data set using different
#' strategies and generates an object of class `mcaps` which stores multiple
#' objects of class [`caps`]. This is a helper function to facilitate comparison
#' of clustering methods and choice of an *optimal* one.
#'
#' @param x A numeric matrix of shape \eqn{N \times M} specifying the grids of
#'   size \eqn{M} on which each of the \eqn{N} curves have been observed.
#' @param y A numeric array of shape \eqn{N \times L \times M} specifying the
#'   \eqn{N}-sample of \eqn{L}-dimensional curves observed on grids of size
#'   \eqn{M}.
#' @param n_clusters_max An integer value specifying the maximum number of
#'   clusers to use. Defaults to `5L`.
#' @param metric A string specifying the metric used to compare curves.
#'   Choices are `"l2"` or `"pearson"`. Defaults to `"l2"`. This is used only
#'   when `warping_class != "srsf"`.
#' @param clustering_method A character vector specifying one or more clustering
#'   methods to be fit. Choices are `"kmeans"`, `"hclust-complete"`,
#'   `"hclust-average"` or `"hclust-single"`. Defaults to all of them.
#' @param warping_class A character vector specifying one or more classes of
#'   warping functions to use for curve alignment. Choices are `"affine"`,
#'   `"dilation"`, `"none"`, `"shift"` or `"srsf"`. Defaults to all of them.
#' @param centroid_type A character vector specifying one or more ways to
#'   compute centroids. Choices are `"mean"` or `"medoid"`. Defaults to all of
#'   them.
#'
#' @return An object of class `mcaps` which is a [`tibble::tibble`] storing the
#'   objects of class [`caps`] in correspondence of each combination of possible
#'   choices from the input arguments.
#' @export
#'
#' @examples
#' out <- compare_caps(
#'   x = simulated30$x,
#'   y = simulated30$y,
#'   n_clusters_max = 2,
#'   warping_class = "affine",
#'   clustering_method = "hclust-complete",
#'   centroid_type = "mean"
#' )
compare_caps <- function(x, y,
                         n_clusters_max = 5L,
                         metric = c("l2", "pearson"),
                         clustering_method = c("kmeans",
                                               "hclust-complete",
                                               "hclust-average",
                                               "hclust-single"),
                         warping_class = c("affine", "dilation", "none",
                                           "shift", "srsf"),
                         centroid_type = c("mean", "medoid", "lowess", "poly")) {
  metric <- rlang::arg_match(metric)
  clustering_method <- rlang::arg_match(clustering_method, multiple = TRUE)
  warping_class <- rlang::arg_match(warping_class, multiple = TRUE)
  centroid_type <- rlang::arg_match(centroid_type, multiple = TRUE)

  df <- tidyr::expand_grid(
    n_clusters = 1:n_clusters_max,
    clustering_method = clustering_method,
    warping_class = warping_class,
    centroid_type = centroid_type
  ) |>
    dplyr::mutate(
      caps_obj = purrr::pmap(
        .l = list(
          .data$n_clusters, .data$clustering_method,
          .data$warping_class, .data$centroid_type
        ),
        .f = \(.n_clusters, .clustering_method, .warping_class, .centroid_type) {
          if (.clustering_method == "kmeans")
            fdakmeans(
              x = x,
              y = y,
              n_clusters = .n_clusters,
              seeding_strategy = "exhaustive-kmeans++",
              warping_class = .warping_class,
              centroid_type = .centroid_type,
              metric = metric,
              use_verbose = FALSE
            )
          else {
            linkage_criterion <- strsplit(.clustering_method, split = "-")[[1]][2]
            fdahclust(
              x = x,
              y = y,
              n_clusters = .n_clusters,
              warping_class = .warping_class,
              centroid_type = .centroid_type,
              metric = metric,
              linkage_criterion = linkage_criterion,
              use_verbose = FALSE
            )
          }
        }
      )
    )

  class(df) <- c("mcaps", class(df))
  df
}

#' Plot for `mcaps` object
#'
#' This is an S3 method implementation of the [`ggplot2::autoplot()`] generic
#' for objects of class `mcaps` to visualize the performances of multiple
#' [`caps`] objects applied on the same data sets either in terms of WSS or in
#' terms of silhouette values.
#'
#' @param object An object of class `mcaps`.
#' @param validation_criterion A string specifying the validation criterion to
#'   be used for the comparison. Choices are `"wss"` or `"silhouette"`. Defaults
#'   to `"wss"`.
#' @param what A string specifying the kind of information to display about the
#'   validation criterion. Choices are `"mean"` (which plots the mean values) or
#'   `"distribution"` (which plots the boxplots). Defaults to `"mean"`.
#' @param ... Other arguments passed to specific methods.
#'
#' @importFrom ggplot2 autoplot
#' @export
#'
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' out <- compare_caps(
#'   x = simulated30$x,
#'   y = simulated30$y,
#'   n_clusters_max = 2,
#'   warping_class = "affine",
#'   clustering_method = "hclust-complete",
#'   centroid_type = "mean"
#' )
#' ggplot2::autoplot(out)
autoplot.mcaps <- function(object,
                           validation_criterion = c("wss", "silhouette"),
                           what = c("mean", "distribution"),
                           ...) {
  validation_criterion <- rlang::arg_match(validation_criterion)
  what <- rlang::arg_match(what)

  if (validation_criterion == "wss") {
    df <- object |>
      dplyr::mutate(
        `Number of clusters` = as.factor(.data$n_clusters),
        value = purrr::map(.data$caps_obj, "distances_to_center")
      ) |>
      dplyr::select(-"caps_obj") |>
      tidyr::unnest(cols = "value") |>
      dplyr::filter(.data$value != 0)
  } else {
    df <- object |>
      dplyr::mutate(
        `Number of clusters` = as.factor(.data$n_clusters),
        value = purrr::map(.data$caps_obj, "silhouettes")
      ) |>
      dplyr::select(-"caps_obj") |>
      tidyr::unnest(cols = "value")
  }

  if (what == "distribution") {
    p <- df |>
      ggplot2::ggplot(ggplot2::aes(
        x = .data$warping_class,
        y = .data$value,
        fill = .data$`Number of clusters`)) +
      ggplot2::geom_boxplot(color = "black") +
      ggplot2::labs(
        title = "Comparison of different clustering strategies",
        subtitle = cli::pluralize("Validation criterion: {toupper(validation_criterion)}"),
        x = "Warping Class",
        y = ""
      )
  } else {
    p <- df |>
      dplyr::group_by(
        .data$`Number of clusters`,
        .data$clustering_method,
        .data$warping_class,
        .data$centroid_type) |>
      dplyr::summarise(value = mean(.data$value)) |>
      dplyr::ungroup() |>
      ggplot2::ggplot(ggplot2::aes(
        x = .data$`Number of clusters`,
        y = .data$value,
        color = .data$warping_class,
        group = .data$warping_class
      )) +
      ggplot2::geom_point() +
      ggplot2::geom_line() +
      ggplot2::labs(
        title = "Comparison of different clustering strategies",
        subtitle = cli::pluralize("Validation criterion: {toupper(validation_criterion)}"),
        y = "",
        color = "Warping Class",
        group = "Warping Class"
      )
  }

  p <- p +
    ggplot2::facet_grid(
      rows = ggplot2::vars(.data$centroid_type),
      cols = ggplot2::vars(.data$clustering_method)
    ) +
    ggplot2::theme_bw()

  if (validation_criterion == "wss")
    p <- p + ggplot2::scale_y_log10()

  p
}

#' Plot for `mcaps` object
#'
#' This is an S3 method implementation of the [`graphics::plot()`] generic for
#' objects of class `mcaps` to visualize the performances of multiple [`caps`]
#' objects applied on the same data sets either in terms of WSS or in terms of
#' silhouette values.
#'
#' @param x An object of class `mcaps`.
#' @inheritParams autoplot.mcaps
#'
#' @return NULL
#'
#' @importFrom graphics plot
#' @export
#'
#' @examples
#' out <- compare_caps(
#'   x = simulated30$x,
#'   y = simulated30$y,
#'   n_clusters_max = 2,
#'   warping_class = "affine",
#'   clustering_method = "hclust-complete",
#'   centroid_type = "mean"
#' )
#' plot(out)
plot.mcaps <- function(x,
                       validation_criterion = c("wss", "silhouette"),
                       what = c("mean", "distribution"),
                       ...) {
  print(autoplot(
    object = x,
    validation_criterion = validation_criterion,
    what = what,
    ...
  ))
}

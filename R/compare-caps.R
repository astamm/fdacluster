#' Generates results of multiple clustering strategies
#'
#' This function searches for clusters in the input data set using different
#' strategies and generates an object of class `mcaps` which stores multiple
#' objects of class [`caps`]. This is a helper function to facilitate comparison
#' of clustering methods and choice of an *optimal* one.
#'
#' @inheritParams fdakmeans
#' @param n_clusters An integer vector specifying a set of clustering partitions
#'   to create. Defaults to `1:5`.
#' @param clustering_method A character vector specifying one or more clustering
#'   methods to be fit. Choices are `"kmeans"`, `"hclust-complete"`,
#'   `"hclust-average"`, `"hclust-single"` or `"dbscan"`. Defaults to all of
#'   them.
#' @param warping_class A character vector specifying one or more classes of
#'   warping functions to use for curve alignment. Choices are `"affine"`,
#'   `"dilation"`, `"none"`, `"shift"` or `"srsf"`. Defaults to all of them.
#' @param centroid_type A character vector specifying one or more ways to
#'   compute centroids. Choices are `"mean"`, `"medoid"`, `"median"`, `"lowess"`
#'   or `"poly"`. Defaults to all of them.
#'
#' @return An object of class `mcaps` which is a [`tibble::tibble`] storing the
#'   objects of class [`caps`] in correspondence of each combination of possible
#'   choices from the input arguments.
#'
#' @export
#' @examples
#' #----------------------------------
#' # Compare k-means results with k = 1, 2, 3, 4, 5 using mean centroid and
#' # various warping classes.
#' \dontrun{
#' sim30_mcaps <- compare_caps(
#'   x = simulated30_sub$x,
#'   y = simulated30_sub$y,
#'   warping_class = c("none", "shift", "dilation", "affine"),
#'   clustering_method = "kmeans",
#'   centroid_type = "mean"
#' )
#' }
#'
#' #----------------------------------
#' # Then visualize the results
#' # Either with ggplot2 via ggplot2::autoplot(sim30_mcaps)
#' # or using graphics::plot()
#' # You can visualize the WSS values:
#' plot(sim30_mcaps, validation_criterion = "wss", what = "mean")
#' plot(sim30_mcaps, validation_criterion = "wss", what = "distribution")
#' # Or the average silhouette values:
#' plot(sim30_mcaps, validation_criterion = "silhouette", what = "mean")
#' plot(sim30_mcaps, validation_criterion = "silhouette", what = "distribution")
compare_caps <- function(x, y,
                         n_clusters = 1:5,
                         metric = c("l2", "pearson"),
                         clustering_method = c("kmeans",
                                               "hclust-complete",
                                               "hclust-average",
                                               "hclust-single",
                                               "dbscan"),
                         warping_class = c("affine", "dilation", "none",
                                           "shift", "srsf"),
                         centroid_type = c("mean", "medoid", "median",
                                           "lowess", "poly"),
                         cluster_on_phase = FALSE) {
  if (!is.numeric(n_clusters))
    cli::cli_abort("The argument {.arg n_clusters} should be an integer/numeric vector.")
  if (length(n_clusters) == 1 && n_clusters == 1)
    cli::cli_abort("It does not make sense to only create a partition with all the data.")
  metric <- rlang::arg_match(metric)
  clustering_method <- rlang::arg_match(clustering_method, multiple = TRUE)
  warping_class <- rlang::arg_match(warping_class, multiple = TRUE)
  centroid_type <- rlang::arg_match(centroid_type, multiple = TRUE)

  df <- tidyr::expand_grid(
    n_clusters = n_clusters,
    clustering_method = clustering_method,
    warping_class = warping_class,
    centroid_type = centroid_type
  ) |>
    dplyr::filter(!(n_clusters > min(n_clusters) & clustering_method == "dbscan")) |>
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
              cluster_on_phase = cluster_on_phase,
              use_verbose = FALSE
            )
          else if (.clustering_method == "dbscan")
            fdadbscan(
              x = x,
              y = y,
              warping_class = .warping_class,
              centroid_type = .centroid_type,
              metric = metric,
              cluster_on_phase = cluster_on_phase,
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
              cluster_on_phase = cluster_on_phase,
              use_verbose = FALSE
            )
          }
        }
      )
    )

  class(df) <- c("mcaps", class(df))
  df
}

#' Visualizes results of multiple clustering strategies using ggplot2
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
#' @return An object of class [`ggplot2::ggplot`].
#'
#' @importFrom ggplot2 autoplot
#' @export
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' p <- ggplot2::autoplot(sim30_mcaps)
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

#' Plots results of multiple clustering strategies
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
#' @examples
#' plot(sim30_mcaps)
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

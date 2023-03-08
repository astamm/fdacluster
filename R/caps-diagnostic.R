#' Diagnostic Plot for [`caps`] Objects
#'
#' This function plots the values of the distance to center and silhouette for
#' each observation. Observations are ordered within cluster by decreasing value
#' of silhouette.
#'
#' @param x An object of class [`caps`].
#'
#' @return An object of class [ggplot2::ggplot].
#' @export
#'
#' @examples
#' out <- fdakmeans(
#'   simulated30$x,
#'   simulated30$y,
#'   seeds = c(1, 21),
#'   n_clusters = 2,
#'   centroid_type = "medoid",
#'   warping_class = "affine",
#'   metric = "pearson"
#' )
#' diagnostic_plot(out)
diagnostic_plot <- function(x) {
  if (!is_caps(x))
    cli::cli_abort("The input argument {.arg x} should be of class {.cls caps}.")

  tibble::tibble(
    `Distance to center` = x$distances_to_center,
    Silhouette = x$silhouettes,
    Membership = as.factor(x$memberships)
  ) |>
    dplyr::mutate(ObservationID = as.character(1:dplyr::n())) |>
    dplyr::group_by(.data$Membership) |>
    dplyr::mutate(ObservationID = forcats::fct_reorder(.data$ObservationID, .data$Silhouette, .desc = TRUE)) |>
    dplyr::ungroup() |>
    tidyr::pivot_longer(cols = c("Distance to center", "Silhouette")) |>
    ggplot2::ggplot(ggplot2::aes(.data$ObservationID, .data$value, fill = .data$Membership)) +
    ggplot2::geom_col(color = "black") +
    ggplot2::facet_wrap(ggplot2::vars(.data$name), ncol = 1, scales = "free") +
    ggplot2::labs(
      title = "CAPS Diagnostic Plot",
      subtitle = cli::pluralize("Number of Clusters: {x$n_clusters} - Metric: {toupper(x$call_args$metric)} - Warping Class: {toupper(x$call_args$warping_class)} - Centroid: {toupper(x$call_args$centroid_type)}"),
      y = ""
    ) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.05))) +
    ggplot2::theme_bw() +
    ggplot2::theme(legend.position = "none")
}

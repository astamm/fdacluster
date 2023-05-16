#' Visualizes the result of a clustering strategy stored in a `caps` object with
#' ggplot2
#'
#' This function creates a visualization of the result of the k-mean alignment
#' algorithm and invisibly returns the corresponding [ggplot2::ggplot] object
#' which enable further customization of the plot. The user can choose to
#' visualize either the amplitude information data in which case original and
#' aligned curves are shown or the phase information data in which case the
#' estimated warping functions are shown.
#'
#' @param object An object of class [`caps`].
#' @param type A string specifying the type of information to display. Choices
#'   are `"amplitude"` for plotting the original and aligned curves which
#'   represent amplitude information data or `"phase"` for plotting the
#'   corresponding warping functions which represent phase information data.
#'   Defaults to `"amplitude"`.
#' @param ... Not used.
#'
#' @return A [ggplot2::ggplot] object invisibly.
#'
#' @export
#' @examplesIf requireNamespace("ggplot2", quietly = TRUE)
#' ggplot2::autoplot(sim30_caps, type = "amplitude")
#' ggplot2::autoplot(sim30_caps, type = "phase")
autoplot.caps <- function(object, type = c("amplitude", "phase"), ...) {
  type <- rlang::arg_match(type)
  if (type == "amplitude") {
    wrangled_data <- plot_data_amplitude(object)
    wrangled_data |>
      ggplot2::ggplot(ggplot2::aes(
        x = .data$grid,
        y = .data$value,
        color = .data$membership,
        group = .data$curve_id
      )) +
      ggplot2::geom_line() +
      ggplot2::facet_grid(
        rows = ggplot2::vars(.data$curve_type),
        cols = ggplot2::vars(.data$component_id)
      ) +
      ggplot2::labs(
        title = "Functional Data",
        subtitle = paste("Class of warping functions:", toupper(object$call_args$warping_class)),
        x = "Observation Grid",
        y = "Component Values"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none")
  } else if (type == "phase") {
    wrangled_data <- plot_data_phase(object)
    wrangled_data |>
      ggplot2::ggplot(ggplot2::aes(
        x = .data$grid,
        y = .data$value,
        color = .data$membership,
        group = .data$curve_id
      )) +
      ggplot2::geom_line() +
      ggplot2::labs(
        title = "Estimated Warping Functions",
        subtitle = paste("Class of warping functions:", toupper(object$call_args$warping_class)),
        x = "Observation Grid",
        y = "Warped Grid Values"
      ) +
      ggplot2::theme_bw() +
      ggplot2::theme(legend.position = "none")
  }
}

#' Plots the result of a clustering strategy stored in a `caps` object
#'
#' This function creates a visualization of the result of the k-mean alignment
#' algorithm **without** returning the plot data as an object. The user can
#' choose to visualize either the amplitude information data in which case
#' original and aligned curves are shown or the phase information data in which
#' case the estimated warping functions are shown.
#'
#' @param x An object of class [`caps`].
#' @inheritParams autoplot.caps
#'
#' @return NULL
#'
#' @export
#' @examples
#' plot(sim30_caps, type = "amplitude")
#' plot(sim30_caps, type = "phase")
plot.caps <- function(x, type = c("amplitude", "phase"), ...) {
  print(autoplot(x, type = type, ...))
}

plot_data_amplitude <- function(x) {
  dplyr::bind_rows(
    wrangle_single_group(x$original_curves, x$grid, x$memberships) |>
      dplyr::mutate(curve_type = "Original Curves"),
    wrangle_single_group(x$aligned_curves, x$grid, x$memberships) |>
      dplyr::mutate(curve_type = "Aligned Curves")
  ) |>
    dplyr::mutate(curve_type = factor(
      .data$curve_type,
      levels = c("Original Curves", "Aligned Curves")
    ))
}

plot_data_phase <- function(x) {
  x$warpings |>
    `colnames<-`(paste0("P", 1:ncol(x$warpings))) |>
    tibble::as_tibble() |>
    dplyr::mutate(
      curve_id = as.factor(1:dplyr::n()),
      membership = as.factor(x$memberships)
    ) |>
    tidyr::pivot_longer(cols = -c("curve_id", "membership")) |>
    tidyr::nest(data = -c("curve_id", "membership")) |>
    dplyr::mutate(grid = purrr::map(.data$membership, \(m) x$grids[m, ])) |>
    tidyr::unnest(cols = c("data", "grid")) |>
    dplyr::select(-"name")
}

wrangle_single_group <- function(curves, grid, memberships) {
  L <- dim(curves)[2]
  1:L |>
    purrr::map(
      .f = wrangle_single_component,
      curves = curves,
      grid = grid,
      memberships = memberships
    ) |>
    purrr::imap(~ dplyr::mutate(.x, component_id = paste("Dimension", .y))) |>
    dplyr::bind_rows()
}

wrangle_single_component <- function(id, curves, grid, memberships) {
  curves <- curves[, id, ]
  N <- nrow(curves)
  M <- ncol(curves)
  K <- length(unique(memberships[memberships > 0]))
  purrr::map(1:K, \(cluster_id) {
    curves[memberships == cluster_id, , drop = FALSE] |>
      `colnames<-`(paste0("P", 1:M)) |>
      tibble::as_tibble() |>
      tidyr::pivot_longer(cols = dplyr::everything()) |>
      tidyr::nest(data = -"name") |>
      dplyr::mutate(grid = grid[cluster_id, ]) |>
      tidyr::unnest(cols = "data") |>
      dplyr::mutate(
        curve_id = which(memberships == cluster_id),
        membership = cluster_id,
        .by = "name") |>
      dplyr::select(-"name")
  }) |>
    dplyr::bind_rows() |>
    dplyr::mutate(
      curve_id = as.factor(.data$curve_id),
      membership = as.factor(.data$membership)
    )
}

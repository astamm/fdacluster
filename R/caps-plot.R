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
    wrangled_data <- plot_amplitude_data(object)
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
    wrangled_data <- plot_phase_data(object)
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

plot_amplitude_data <- function(x) {
  original_data <- format_viz(x$original_grids, x$original_curves, x$memberships)
  aligned_data <- format_viz(x$aligned_grids, x$original_curves, x$memberships)
  original_data$curve_type <- "Original Curves"
  aligned_data$curve_type <- "Aligned Curves"
  curve_type_levels <- c("Original Curves", "Aligned Curves")
  original_data$curve_type <- factor(original_data$curve_type, levels = curve_type_levels)
  aligned_data$curve_type <- factor(aligned_data$curve_type, levels = curve_type_levels)
  rbind(original_data, aligned_data)
}

plot_phase_data <- function(x) {
  data.frame(
    grid = as.vector(x$original_grids),
    value = as.vector(x$aligned_grids),
    curve_id = as.factor(1:nrow(x$original_grids)),
    membership = as.factor(x$memberships)
  )
}

format_viz <- function(grids, curves, memberships) {
  dims <- dim(curves)
  N <- dims[1]
  L <- dims[2]
  M <- dims[3]
  out <- lapply(1:L, function(l) {
    unicurves <- curves[, l, ]
    data.frame(
      grid = as.vector(grids),
      value = as.vector(unicurves),
      membership = memberships,
      curve_id = 1:N,
      component_id = l
    )
  }) |>
    do.call(what = rbind, args = _)
  out$membership <- as.factor(out$membership)
  out$curve_id <- as.factor(out$curve_id)
  out$component_id <- as.factor(out$component_id)
  out
}

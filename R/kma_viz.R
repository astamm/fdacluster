#' Plot for \code{kmap} objects
#'
#' @param x The [kma] object to be plotted.
#' @param type A string specifying the type of information to display. Choices
#'   are \code{"data"} for plotting the original and aligned curves (default) or
#'   \code{"warping"} for plotting the corresponding warping functions.
#' @param number_of_displayed_points The number of points to used for display.
#'   It is set as the minimum between this parameter and the number of points in
#'   the original data set. Defaults to 50.
#' @param ... Other graphical parameters (see [par][graphics::par]). Ignored for
#'   now.
#'
#' @return A [ggplot][ggplot2::ggplot] object invisibly.
#' @export
#'
#' @examples
#' res <- kma(
#'   simulated30$x,
#'   simulated30$y,
#'   seeds = c(1, 21),
#'   n_clust = 2,
#'   center_method = "medoid",
#'   warping_method = "affine",
#'   dissimilarity_method = "pearson"
#' )
#'
#' plot(res, type = "data")
#' plot(res, type = "warping")
plot.kma <- function(x, type = "data", number_of_displayed_points = 50, ...) {
  if (type == "data")
    plot_data(x, number_of_displayed_points)
  else if (type == "warping")
    plot_warping(x, number_of_displayed_points)
  else
    stop("Unsupported type of display for kma objects.")
}

plot_data <- function(obj, type = "data", number_of_displayed_points = 50) {
  n <- dim(obj$y)[1]
  d <- dim(obj$y)[2]
  p <- dim(obj$y)[3]

  original_grids <- 1:p %>%
    purrr::map(~ obj$x[, .x, drop = FALSE]) %>%
    purrr::set_names(paste0("P", 1:p)) %>%
    as_tibble() %>%
    dplyr::mutate(curve_id = 1:dplyr::n()) %>%
    tidyr::pivot_longer(
      cols = -.data$curve_id,
      names_to = "point_id",
      values_to = "grid"
    ) %>%
    dplyr::mutate(type = "Original Curves")

  warped_grids <- 1:p %>%
    purrr::map(~ obj$x_final[, .x, drop = FALSE]) %>%
    purrr::set_names(paste0("P", 1:p)) %>%
    as_tibble() %>%
    dplyr::mutate(curve_id = 1:dplyr::n()) %>%
    tidyr::pivot_longer(
      cols = -.data$curve_id,
      names_to = "point_id",
      values_to = "grid"
    ) %>%
    dplyr::mutate(type = "Aligned Curves")

  original_values <- 1:p %>%
    purrr::map(~ {
      df <- matrix(obj$y[, , .x], ncol = d)
      colnames(df) <- paste("Dimension", 1:d)
      as_tibble(df) %>%
        dplyr::mutate(curve_id = 1:dplyr::n())
    }) %>%
    purrr::set_names(paste0("P", 1:p)) %>%
    dplyr::bind_rows(.id = "point_id") %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("Dim"),
      names_to = "dimension_id",
      values_to = "value"
    )

  center_grids <- 1:p %>%
    purrr::map(~ obj$x_centers_final[, .x, drop = FALSE]) %>%
    purrr::set_names(paste0("P", 1:p)) %>%
    as_tibble() %>%
    dplyr::mutate(curve_id = 1:dplyr::n()) %>%
    tidyr::pivot_longer(
      cols = -.data$curve_id,
      names_to = "point_id",
      values_to = "grid"
    ) %>%
    dplyr::mutate(type = "Aligned Curves")

  center_values <- 1:p %>%
    purrr::map(~ {
      df <- matrix(obj$y_centers_final[, , .x], ncol = d)
      colnames(df) <- paste("Dimension", 1:d)
      as_tibble(df) %>%
        dplyr::mutate(curve_id = 1:dplyr::n())
    }) %>%
    purrr::set_names(paste0("P", 1:p)) %>%
    dplyr::bind_rows(.id = "point_id") %>%
    tidyr::pivot_longer(
      cols = dplyr::starts_with("Dim"),
      names_to = "dimension_id",
      values_to = "value"
    )

  number_of_displayed_points <- min(number_of_displayed_points, p)

  df <- dplyr::bind_rows(
    original_values %>%
      dplyr::left_join(original_grids, by = c("point_id", "curve_id")),
    original_values %>%
      dplyr::left_join(warped_grids, by = c("point_id", "curve_id"))
  ) %>%
    dplyr::mutate(type = factor(.data$type, c("Original Curves", "Aligned Curves"))) %>%
    dplyr::left_join(tibble(
      curve_id = 1:n,
      membership = as.factor(obj$labels)
    ), by = "curve_id") %>%
    dplyr::group_by(.data$curve_id, .data$dimension_id, .data$type) %>%
    dplyr::slice(seq(
      from = 1,
      to = dplyr::n(),
      by = round((dplyr::n() - 1) / number_of_displayed_points)
    )) %>%
    dplyr::ungroup()

  df_mean <- center_values %>%
    dplyr::left_join(center_grids, by = c("point_id", "curve_id")) %>%
    dplyr::mutate(
      membership = as.factor(.data$curve_id),
      type = factor(.data$type, c("Original Curves", "Aligned Curves"))
    ) %>%
    dplyr::group_by(.data$curve_id, .data$dimension_id, .data$type) %>%
    dplyr::slice(seq(
      from = 1,
      to = dplyr::n(),
      by = round((dplyr::n() - 1) / number_of_displayed_points)
    )) %>%
    dplyr::ungroup()

  df %>%
    ggplot(aes(.data$grid, .data$value, group = .data$curve_id, color = .data$membership)) +
    geom_line(alpha = 0.3) +
    geom_line(data = df_mean, size = 1.5) +
    facet_wrap(vars(.data$type, .data$dimension_id), nrow = 2, scales = "free") +
    theme_bw() +
    theme(legend.position = "top") +
    labs(
      title = "Functional Data",
      subtitle = paste("Class of warping functions:", toupper(obj$warping_method)),
      x = "Grid",
      y = "Values",
      color = "Group membership"
    )
}

plot_warping <- function(obj, number_of_displayed_points = 50) {
  if (obj$warping_method == "affine")
    df <- data_affine(obj, number_of_displayed_points)
  else if (obj$warping_method == "dilation")
    df <- data_dilation(obj, number_of_displayed_points)
  else if (obj$warping_method == "shift")
    df <- data_shift(obj, number_of_displayed_points)
  else
    stop("Unsupported warping family for display.")
  df %>%
    tidyr::unnest(cols = .data$x:.data$y) %>%
    ggplot(aes(.data$x, .data$y, color = .data$membership, group = .data$id)) +
    geom_line() +
    theme_bw() +
    theme(legend.position = "top") +
    labs(
      title = "Estimated Warping Functions",
      subtitle = paste("Class of warping functions:", toupper(obj$warping_method)),
      x = "Original grids",
      y = "Warped grids",
      color = "Group membership"
    )
}

data_affine <- function(obj, number_of_displayed_points = 50) {
  obj$parameters %>%
    `colnames<-`(c("slope", "intercept")) %>%
    as_tibble() %>%
    dplyr::mutate(
      x = purrr::map2(
        .x = .data$slope,
        .y = .data$intercept,
        .f = ~ seq(0, 1, length.out = number_of_displayed_points)
      ),
      y = purrr::map2(
        .x = .data$slope,
        .y = .data$intercept,
        .f = ~ .x * seq(0, 1, length.out = number_of_displayed_points) + .y
      ),
      id = 1:dplyr::n(),
      membership = as.factor(obj$labels)
    )
}

data_dilation <- function(obj, number_of_displayed_points = 50) {
  obj$parameters %>%
    `colnames<-`("slope") %>%
    as_tibble() %>%
    dplyr::mutate(
      x = purrr::map(
        .x = .data$slope,
        .f = ~ seq(0, 1, length.out = number_of_displayed_points)
      ),
      y = purrr::map(
        .x = .data$slope,
        .f = ~ .x * seq(0, 1, length.out = number_of_displayed_points)
      ),
      id = 1:dplyr::n(),
      membership = as.factor(obj$labels)
    )
}

data_shift <- function(obj, number_of_displayed_points = 50) {
  obj$parameters %>%
    `colnames<-`("intercept") %>%
    as_tibble() %>%
    dplyr::mutate(
      x = purrr::map(
        .x = .data$intercept,
        .f = ~ seq(0, 1, length.out = number_of_displayed_points)
      ),
      y = purrr::map(
        .x = .data$intercept,
        .f = ~ seq(0, 1, length.out = number_of_displayed_points) + .x
      ),
      id = 1:dplyr::n(),
      membership = as.factor(obj$labels)
    )
}

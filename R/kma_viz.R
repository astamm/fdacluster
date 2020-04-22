#' Plot for \code{kmap} objects
#'
#' @param x The \code{kma} object to be plotted.
#' @param number_of_displayed_points The number of points to graphically
#'   represents. It is set as the minimum between this parameter and the number
#'   of points in the original data set. Defaults to 50.
#' @param ... Additional parameters passed to the \code{\link[graphics]{plot}}
#'   function.
#'
#' @return A `ggplot` object invisibly.
#' @export
#'
#' @examples
#' res <- kma(
#'   simulated30$x,
#'   simulated30$y,
#'   seeds = c(1, 11, 21),
#'   n_clust = 3,
#'   center_method = "medoid",
#'   warping_method = "affine",
#'   dissimilarity_method = "pearson"
#' )
#'
#' plot(res)
plot.kma <- function(x, number_of_displayed_points = 50, ...) {
  n <- dim(x$y)[1]
  d <- dim(x$y)[2]
  p <- dim(x$y)[3]

  original_grids <- 1:p %>%
    purrr::map(~ x$x[, .x, drop=FALSE]) %>%
    purrr::set_names(paste0("P", 1:p)) %>%
    as_tibble() %>%
    dplyr::mutate(curve_id = 1:dplyr::n()) %>%
    tidyr::pivot_longer(
      cols = -curve_id,
      names_to = "point_id",
      values_to = "grid"
    ) %>%
    dplyr::mutate(type = "Original Curves")

  warped_grids <- 1:p %>%
    purrr::map(~ x$x_final[, .x, drop=FALSE]) %>%
    purrr::set_names(paste0("P", 1:p)) %>%
    as_tibble() %>%
    dplyr::mutate(curve_id = 1:dplyr::n()) %>%
    tidyr::pivot_longer(
      cols = -curve_id,
      names_to = "point_id",
      values_to = "grid"
    ) %>%
    dplyr::mutate(type = "Aligned Curves")

  original_values <- 1:p %>%
    purrr::map(~ {
      df <- as.matrix(x$y[, , .x, drop=FALSE])
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
    purrr::map(~ x$x_centers_final[, .x, drop=FALSE]) %>%
    purrr::set_names(paste0("P", 1:p)) %>%
    as_tibble() %>%
    dplyr::mutate(curve_id = 1:dplyr::n()) %>%
    tidyr::pivot_longer(
      cols = -curve_id,
      names_to = "point_id",
      values_to = "grid"
    ) %>%
    dplyr::mutate(type = "Aligned Curves")

  center_values <- 1:p %>%
    purrr::map(~ {
      df <- as.matrix(x$y_centers_final[, , .x, drop=FALSE])
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
    dplyr::mutate(type = factor(type, c("Original Curves", "Aligned Curves"))) %>%
    dplyr::left_join(tibble(
      curve_id = 1:n,
      membership = as.factor(x$labels)
    ), by = "curve_id") %>%
    dplyr::group_by(curve_id, dimension_id, type) %>%
    dplyr::slice(seq(1, dplyr::n(), by = (dplyr::n() - 1) / number_of_displayed_points)) %>%
    dplyr::ungroup()

  df_mean <- center_values %>%
    dplyr::left_join(center_grids, by = c("point_id", "curve_id")) %>%
    dplyr::mutate(
      membership = as.factor(curve_id),
      type = factor(type, c("Original Curves", "Aligned Curves"))
    ) %>%
    dplyr::group_by(curve_id, dimension_id, type) %>%
    dplyr::slice(seq(1, dplyr::n(), by = (dplyr::n() - 1) / number_of_displayed_points)) %>%
    dplyr::ungroup()

  df %>%
    ggplot(aes(grid, value, group = curve_id, color = membership)) +
    geom_line(alpha = 0.3) +
    geom_line(data = df_mean, size = 1.5) +
    facet_grid(rows = vars(type), cols = vars(dimension_id)) +
    theme_bw() +
    theme(legend.position = "top")
}

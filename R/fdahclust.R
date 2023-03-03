#' HAC
#'
#' blabla
#'
#' @inheritParams fdakmeans
#'
#' @return An object of class [`caps`].
#'
#' @export
#' @examples
#' out <- fdahclust(simulated30$x, simulated30$y)
fdahclust <- function(x, y,
                      n_clusters = 1L,
                      warping_class = c("affine", "dilation", "none", "shift", "srsf"),
                      seeds = NULL,
                      maximum_number_of_iterations = 100L,
                      centroid_type = c("mean", "medoid"),
                      metric = c("l2", "pearson"),
                      warping_options = c(0.15, 0.15),
                      number_of_threads = 1L,
                      parallel_method = 0L,
                      distance_relative_tolerance = 0.001,
                      use_fence = FALSE,
                      check_total_dissimilarity = TRUE,
                      use_verbose = TRUE,
                      compute_overall_center = FALSE) {

}

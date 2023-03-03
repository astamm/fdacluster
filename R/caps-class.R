#' Class for Clustering with Amplitude and Phase Separation
#'
#' The k-means algorithm with joint amplitude and phase separation produces a
#' number of outputs. This class is meant to encapsulate them into a single
#' object for providing dedicated `S3` methods for e.g. plotting, summarizing,
#' etc. The name of the class stems from **C**lustering with **A**mplitude and
#' **P**hase **S**eparation.
#'
#' An object of class [`caps`] is a list with the following components:
#'
#' - `original_curves`: A numeric matrix of shape \eqn{N \times L \times M}
#' storing a sample with the \eqn{N} \eqn{L}-dimensional original curves
#' observed on grids of size \eqn{M}.
#' - `aligned_grids`: A numeric matrix of shape \eqn{N \times L \times M}
#' storing a sample with the \eqn{N} \eqn{L}-dimensional aligned curves observed
#' on grids of size \eqn{M}.
#' - `center_curves`: A numeric matrix of shape \eqn{K \times L \times M}
#' storing the \eqn{K} centers which are \eqn{L}-dimensional curves observed on
#' a grid of size \eqn{M};
#' - `grid`: A numeric vector of length \eqn{M} storing the
#' common grid of size \eqn{M} on which curves have been observed;
#' - `n_clusters`: An integer value storing the number of clusters;
#' - `memberships`: An integer vector of length \eqn{N} storing the cluster ID
#' which each curve belongs to;
#' - `distances_to_center`: A numeric vector of length \eqn{N} storing the
#' distance of each curve to the center of its cluster;
#' - `warpings`: A numeric matrix of shape \eqn{N \times M} storing the
#' estimated warping functions for each of the \eqn{N} curves evaluted on the
#' common `grid` of size \eqn{M};
#' - `n_iterations`: An integer value storing the number of iterations
#' performed until convergence;
#' - `call_name`: A string storing the name of the function that was used to
#' produce the k-means alignment results;
#' - `call_args`: A list containing the exact arguments that were passed to
#' the function `call_name` that produced this output.
#'
#' @param x A list coercible into an object of class [`caps`].
#'
#' @name caps
#'
#' @examples
#' res <- fdakmeans(
#'   simulated30$x,
#'   simulated30$y,
#'   seeds = c(1, 21),
#'   n_clusters = 2,
#'   centroid_type = "medoid",
#'   warping_class = "affine",
#'   metric = "pearson"
#' )
#' as_caps(res)
#' is_caps(res)
NULL

#' @rdname caps
#' @export
as_caps <- function(x) {
  if (is_caps(x)) return(x)

  if (!inherits(x, "list"))
    cli::cli_abort("The input argument {.arg x} should be a list.")

  if (length(x) != 11)
    cli::cli_abort("The input argument {.arg x} should be a list of length 11.")

  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  if (any(names(x) != expected_names))
    cli::cli_abort("The input argument {.arg x} should be a list with components {expected_names}.")

  class(x) <- c("caps", class(x))

  x
}

#' @rdname caps
#' @export
is_caps <- function(x) {
  "caps" %in% class(x)
}

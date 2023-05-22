#' Class for clustering with amplitude and phase separation
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
#' - `aligned_curves`: A numeric matrix of shape \eqn{N \times L \times M}
#' storing a sample with the \eqn{N} \eqn{L}-dimensional aligned curves observed
#' on grids of size \eqn{M}.
#' - `center_curves`: A numeric matrix of shape \eqn{K \times L \times M}
#' storing the \eqn{K} centers which are \eqn{L}-dimensional curves observed on
#' a grid of size \eqn{M};
#' - `warpings`: A numeric matrix of shape \eqn{N \times M} storing the
#' estimated warping functions for each of the \eqn{N} curves evaluated on the
#' within-cluster common `grids` of size \eqn{M};
#' - `grids`: A numeric matrix of size \eqn{K \times M} storing the
#' within-cluster common grids of size \eqn{M} on which curves are evaluated;
#' - `n_clusters`: An integer value storing the number of clusters;
#' - `memberships`: An integer vector of length \eqn{N} storing the cluster ID
#' which each curve belongs to;
#' - `distances_to_center`: A numeric vector of length \eqn{N} storing the
#' distance of each curve to the center of its cluster;
#' - `silhouettes`: A numeric vector of length \eqn{N} storing the silhouette
#' values of each observation;
#' - `n_iterations`: An integer value storing the number of iterations
#' performed until convergence;
#' - `call_name`: A string storing the name of the function that was used to
#' produce the k-means alignment results;
#' - `call_args`: A list containing the exact arguments that were passed to
#' the function `call_name` that produced this output.
#'
#' @param x A list coercible into an object of class [`caps`].
#'
#' @return The function [`as_caps()`] returns an object of class [`caps`]. The
#'   function [`is_caps()`] returns a boolean which evaluates to `TRUE` is the
#'   input object is of class [`caps`].
#'
#' @name caps
#'
#' @examples
#' as_caps(sim30_caps)
#' is_caps(sim30_caps)
NULL

#' @rdname caps
#' @export
as_caps <- function(x) {
  if (is_caps(x)) return(x)

  if (!inherits(x, "list"))
    cli::cli_abort("The input argument {.arg x} should be a list.")

  if (length(x) != 14)
    cli::cli_abort("The input argument {.arg x} should be a list of length 14.")

  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
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

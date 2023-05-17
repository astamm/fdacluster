#' Simulated data for examples
#'
#' A data set containing 30 simulated uni-dimensional curves.
#'
#' @format A list with abscissas x and values y:
#' \describe{
#'   \item{x}{Matrix 30x200;}
#'   \item{y}{Array 30x1x200.}
#' }
#'
"simulated30"

#' Simulated data for examples
#'
#' A data set containing 30 simulated uni-dimensional curves.
#'
#' @format A list with abscissas x and values y:
#' \describe{
#'   \item{x}{Matrix 30x30;}
#'   \item{y}{Array 30x1x30.}
#' }
#'
"simulated30_sub"

#' A `caps` object from simulated data for examples
#'
#' An object of class [`caps`] storing the result of the [`fdakmeans()`]
#' function applied on the data set [`simulated30`] using the affine warping
#' class and the Pearson metric and searching for 2 clusters.
#'
#' @format An object of class [`caps`].
#'
"sim30_caps"

#' An `mcaps` object from simulated data for examples
#'
#' An object of class `mcaps` storing the result of the [`compare_caps()`]
#' function applied on the data set [`simulated30_sub`] for comparing the
#' clustering structures found by the [`fdakmeans()`] function with `mean`
#' centroid type used with various classes of warping functions and varying
#' number of clusters.
#'
#' @format An object of class `mcaps` which is effectively a [tibble::tibble]
#'   with 5 columns and as many rows as there are clustering strategies to
#'   compare. The 5 column-variables are:
#'
#'   - `n_clusters`: The number of clusters;
#'   - `clustering_method`: The clustering method;
#'   - `warping_class`: The class of warping functions used for curve alignment;
#'   - `centroid_type`: The type of centroid used to compute a cluster
#'   representative;
#'   - `caps_obj`: The result of the corresponding clustering strategy as
#'   objects of class [`caps`].
#'
"sim30_mcaps"

#' Simulated data from the CSDA paper
#'
#' A data set containing 90 simulated uni-dimensional curves.
#'
#' @format A list with abscissas x and values y:
#' \describe{
#'   \item{x}{Vector of size 100;}
#'   \item{y}{Matrix if size 90x100.}
#' }
#'
"simulated90"

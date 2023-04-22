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
#'   \item{x}{Matrix 30x50;}
#'   \item{y}{Array 30x1x50.}
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
#' function applied on the data set [`simulated30`] for comparing the clustering
#' structures found by the [`fdakmeans()`] function with `mean` centroid type
#' used with various classes of warping functions and varying number of
#' clusters.
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

#' Subset of the AneuRisk65 benchmark data set
#'
#' A data set containing the first derivative of the three-dimensional
#' coordinates of the centerline of the internal carotid artery of 65 patients.
#'
#' @format A list with 2 components:
#' \describe{
#'   \item{x}{A 65 x 1380 matrix containing, in each row, the evaluation grid
#'   for each patient;}
#'   \item{y}{A 65 x 3 x 1380 array containing, in each row, the values of the
#'   first derivative of each of the 3D coordinates of the ICA centerline,
#'   stored by a row in a matrix.}
#' }
#'
#' @source This is a subset of the [AneuRisk65 benchmark data
#'   set](https://mox.polimi.it/research-areas/statistics/) provided by the AneuRisk
#'   project.
#'
#' @references \enumerate{
#'   \item Sangalli, L.M., Secchi, P. and Vantini, S. (2014), [AneuRisk65: A
#'   dataset of three-dimensional cerebral vascular
#'   geometries](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-8/issue-2/AneuRisk65-A-dataset-of-three-dimensional-cerebral-vascular-geometries/10.1214/14-EJS938.full),
#'   Electronic Journal of Statistics, 8 (2), 1879-1890.
#'   \item Sangalli, L.M., Secchi, P. and Vantini, S. (2014), [Analysis of
#'   AneuRisk65 data: K-mean
#'   Alignment](https://projecteuclid.org/journals/electronic-journal-of-statistics/volume-8/issue-2/Analysis-of-AneuRisk65-data-k-mean-alignment/10.1214/14-EJS938A.full),
#'   Electronic Journal of Statistics, 8 (2), 1891-1904.
#' }
#'
"aneurisk65"

#' @import nloptr
#' @importFrom ggplot2 autoplot
#' @importFrom graphics plot
#' @importFrom Rcpp sourceCpp
#' @importFrom rlang .data
#' @useDynLib fdacluster, .registration = TRUE
#'
#' @references Arthur, D., and S. Vassilvitskii. 2007. “K-Means++ the Advantages
#'   of Careful Seeding.” In Proceedings of the Eighteenth Annual ACM-SIAM
#'   Symposium on Discrete Algorithms, 1027–35.
#' @references Marron, J. S., J. O. Ramsay, L. M. Sangalli, and A. Srivastava.
#'   2014. “Statistics of Time Warpings and Phase Variations.”
#' @references Marron, J. S., J. O. Ramsay, L. M. Sangalli, and A. Srivastava.
#'   2015. “Functional Data Analysis of Amplitude and Phase Variation.”
#'   Statistical Science, 468–84.
#' @references Ramsay, J., and B. W. Silverman. 2005. Functional Data Analysis.
#'   Springer Series in Statistics. Springer.
#' @references Sangalli, L. M., P. Secchi, S. Vantini, and V. Vitelli. 2010.
#'   “K-Mean Alignment for Curve Clustering.” Computational Statistics & Data
#'   Analysis 54 (5): 1219–33.
#' @references Tucker, J. D., W. Wu, and A. Srivastava. 2013. “Generative Models
#'   for Functional Data Using Phase and Amplitude Separation.” Computational
#'   Statistics & Data Analysis 61: 50–66.
#' @references Vantini, S. 2012. “On the Definition of Phase and Amplitude
#'   Variability in Functional Data Analysis.” Test 21 (4): 676–96.
#'   https://doi.org/10.1007/s11749-011-0268-9.
#'
#' @keywords internal
"_PACKAGE"

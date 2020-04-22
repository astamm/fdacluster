#' K-mean alignment algorithm and variants for functional data
#'
#' The **fdakmapp** package allows to jointly perform clustering
#' and alignment of functional data.
#'
#' @references
#' \enumerate{
#'   \item Sangalli, L.M., Secchi, P., Vantini, S. and Vitelli, V. (2010),
#'   [K-mean alignment for curve
#'   clustering](https://www.sciencedirect.com/science/article/abs/pii/S0167947309004605),
#'   Computational Statistics and Data Analysis, 54, 1219-1233.
#'   \item Sangalli, L.M., Secchi, P. and Vantini, S. (2014), [Analysis of
#'   AneuRisk65 data: K-mean
#'   Alignment](http://projecteuclid.org/euclid.ejs/1414588177), Electronic
#'   Journal of Statistics, 8 (2), 1891-1904.
#' }
#'
#' @name fdakmapp
#' @docType package
#' @useDynLib fdakmapp, .registration=TRUE
#' @import nloptr
#' @importFrom Rcpp sourceCpp
NULL

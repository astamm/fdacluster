#'@title Functional Data Analysis Plus: K-Mean/Medoid Alignment:
#'@description Fdakmapp is a package that allows to jointly perform clustering
#' and alignment of a functional dataset (multidimensional or unidimensional functions).
#'
#'@references
#'\itemize{
#'\item Sangalli, L.M., Secchi, P., Vantini, S., Vitelli, V., 2010. \emph{"K-mean alignment for curve clustering"}. Computational Statistics and Data Analysis, 54, 1219-1233.
#'\item Sangalli, L.M., Secchi, P., Vantini, S., 2014. \emph{"Analysis of AneuRisk65 data: K-mean Alignment"}. Electronic Journal of Statistics, Special Section on "Statistics of Time Warpings and Phase Variations", Vol. 8, No. 2, 1891-1904.
#'}
#'@seealso
#'\code{\link{kmap}}.
#'
#'@examples
#' ############# EXE ##########################
#' res<-kmap(x=aneurisk65$x, y=aneurisk65$y, n_clust=2,seeds=c(32,64))
#'
#' ############# OUTPUT################
#' kmap_show_results(res,TRUE,TRUE)
#'
#'@docType package
#'@name fdakmapp
NULL





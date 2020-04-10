#' K-mean alignment and variants for functional data
#'
#' This function jointly performs clustering and alignment of functional data
#' sets with possibly multidimensional domains.
#'
#' @param x Numeric matrix [\emph{nObs} x \emph{nPts}] or vector [\emph{nPts}]:
#'   the grid on which where each function is evaluated, where \emph{nObs} is
#'   the number of observations and \emph{nPts} is the grid size. Grids are
#'   allowed to be different for each observation. If \code{x} is a vector of
#'   length \emph{nPts}, it is assumed to be the common grid on which all
#'   observations are evaluated. Smaller grids are completed with missing values
#'   to match the size of the largest.
#' @param y Numeric matrix [\emph{nObs} x \emph{nPts}] or array [\emph{nObs} x
#'   \emph{nPts} x \emph{nDim}]: values of the observed functions on the grid
#'   \code{x}, where \emph{d} is the dimension of the domain of definition of
#'   the functions.
#' @param seeds Numeric vector [nClusters]: indices of the observations used as
#'   initial centers for the k-mean algorithm. If \code{seeds = NULL} (default),
#'   they are randomly chosen among the \emph{nObs} original functions.
#' @param n_clust Scalar: number of clusters to look for (default: 1).
#' @param warping_method Character: class of warping functions among
#'   \code{noalign} for no alignment which is equivalent to performing k-means
#'   clustering only, \code{shift}, \code{dilation} and \code{affine} (default).
#' @param center_method Character: method to find the cluster representative
#'   among \code{mean} (default), \code{medoid}, \code{pseudomedoid} and
#'   \code{median}.
#' @param similarity_method Character: metric to measure dissimilarity between
#'   two observations. Choices are \code{pearson} (default), \code{l2},
#'   \code{l2weighted} and \code{l2first}.
#' @param optim_method Character: optimization method for searching the optimal
#'   warping functions at each iteration. The only available method at the
#'   moment is \code{bobyqa}.
#' @param warping_opt Numeric vector. It depends on the class of warping
#'   functions. For affine functions, it is a 2D vector specifying the upper
#'   bounds for dilation and shift. For dilations, a scalar specifying the upper
#'   bound for dilation. For shifts, a scalar specifying the upper bound for
#'   shift. For no alignment, it is an empty vector. Default value is
#'   \code{c(0.15, 0.15)} for affine warping functions.
#' @param center_opt Numeric vector. It depends on the cluster representative
#'   method. It is an empty vector for all methods except \code{mean}. In this
#'   latter case, it is a 2D vector of the form \code{c(span, delta)} specifying
#'   the parameters for the \code{lowess} approximation. Default value is
#'   \code{c(0.01, 0.1)}.
#' @param out_opt Numeric vector. It is a 3D vector of the form \code{c(nOut,
#'   tol, maxIter)}. It contains the size of the output grid on which the
#'   centers will be evaluated and two stopping criteria for the optimization
#'   algorithm. The defaut value is \code{c(100, 0.001, 100)}.
#' @param fence Boolean. When \code{TRUE}, a control is activated at the end of
#'   each iteration. The aim of the control is to avoid warping outliers with
#'   respect to their computed distributions. Default value is \code{FALSE}
#'   because it is time-consuming.
#' @param check_total_similarity Boolean. When \code{TRUE} (default), at each
#'   iteration, the algorithm checks if the total similarity is improving and
#'   stops if it is not.
#' @param show_iter Boolean. When \code{TRUE} (default), information pertaining
#'   to the current iteration is displayed in the console.
#' @param comp_original_center Boolean. When \code{TRUE}, the initial center
#'   with relative dissimilarities is computed. Otherwise, this step is skipped.
#'   Default value is \code{FALSE} because it is time-consuming.
#' @param par_opt Numeric vector. Parallel options as a 2D vector of the form
#'   \code{c(nThreads, parVersion)}, where \code{nThreads} specifies the number
#'   of threads to run on parallel and \code{parVersion} is either 0 for a
#'   trivial parallelization in which each thread computes the center of a
#'   cluster or 1 for a more efficient parallelization in which all threads
#'   compute the centers sequentially. The latter parallelizatiobn mode is
#'   available only when \code{center_method = 'medoid'}. Default value is
#'   \code{c(1, 0)}.
#'
#' @return The function output is a \code{kmap} object, which is a list with the following elements:
#' \item{x}{ as input.}
#' \item{y}{ as input. }
#' \item{seeds}{ vector with the indeces used in the algorithm.}
#' \item{warping_method}{ as input.}
#' \item{similarity_method}{ as input. }
#' \item{center_method}{ as input. }
#' \item{iterations}{scalar: total number of iterations performed by kma function.}
#' \item{n_clust}{ as input. }
#' \item{x.center.orig }{numeric vector \emph{n_out}: abscissa of the center computed if \emph{comp_original_center}=TRUE.}
#' \item{y.center.orig }{numeric vector \emph{n_out} or matrix \emph{n_out} X \emph{n_dim}: value of the center computed if \emph{comp_original_center}=TRUE.}
#' \item{similarity.origin}{numeric vector \emph{n_obs} dissimilarity,similarity or distance of the original center respect the obserbations computed if \emph{comp_original_center}=TRUE.}
#' \item{x.final }{matrix [n.func X grid.size]: aligned abscissas.}
#' \item{n.clust.final }{ scalar: final number of clusters. Note that n.clust.final may differ from initial number of clusters (i.e.,from n.clust) if some clusters are found to be empty.}
#' \item{x.centers.final }{matrix [n.clust.final X grid.size]: abscissas of the final function centers.}
#' \item{y.centers.final }{matrix [n.clust.final X n.out] or array [n.clust.final X n.out x n_dim] , contain the evaluations of the final functions centers.}
#' \item{templates_vec }{list iteration : each element of the list contain centers of that iteration.}
#' \item{x_out_vec }{list iteration : each element of the list contain the abscissa of the centers of that iteration.}
#' \item{labels}{vector n_obs: cluster assignments.}
#' \item{similarity.final}{vector [n_obs]: similarities,dissimilarities or distance between each function and the center of the cluster the function is assigned to.}
#' \item{parameters.list}{list [iterations]: warping parameters at each iteration.}
#' \item{parameters}{matrix [n_par X n_obs]: warping parameters applied to the original abscissas x to obtain the aligned abscissas x.final.}
#' \item{timer}{vector: time of execution by step. }
#' @export
#'
#' @examples
kmap <- function(x, y,
                 seeds = NULL,
                 n_clust = 1,
                 warping_method = 'affine',
                 center_method = 'mean',
                 similarity_method = 'pearson',
                 optim_method = 'bobyqa',
                 warping_opt = c(0.15, 0.15),
                 center_opt = c(0.01, 0.1),
                 out_opt = c(100, 0.001, 100),
                 fence = FALSE,
                 check_total_similarity = TRUE,
                 show_iter = TRUE,
                 comp_original_center = FALSE,
                 par_opt = c(1, 0))
{

  if (is.null(y))
    stop("Provide a valid function")

  if (length(dim(y)) == 2) {
    y <- array(y, c(dim(y)[1], dim(y)[2], 1))
  }

  if (is.vector(x))
   x <- t(replicate(dim(y)[1],x,1))

  if(is.null(seeds)){
    nseeds<-sample(0:(nrow(y)-1),n_clust)
  }else{
    nseeds<-(seeds-1)
  }
  if ((sum(nseeds<0)+sum(nseeds>=nrow(y)))!=0 )
    stop("seeds indexes have to be in observations range")

  out<-.Call('_fdakmapp_kmap', PACKAGE = 'fdakmapp', x, y, nseeds, n_clust, warping_method, center_method, similarity_method, optim_method, warping_opt, center_opt, out_opt, fence, check_total_similarity, show_iter, comp_original_center,par_opt)

  ## gestione timer  ################################################
  time<-diff(round(out$timer/1000000000,4))
  t<-data.frame(0,0,0,0,0,0)
  names(t)<-c("start","warping","fece/norm","templates","output","total")
  rownames(t)<-c("sec")
  t[1]<-time[1]
  for(i in 0:(out$iterations-1)){
    t[2]=t[2]+time[2+(i*3)]
    t[3]=t[3]+time[3+(i*3)]
    t[4]=t[4]+time[4+(i*3)]
  }
  t[5]=time[out$iterations*3+2]
  t[6]=out$timer[length(out$timer)]/1000000000
  out$timer<-t
    #######################################################################

  out <- c(
    x = list(x),
    y = list(y),
    seeds = list(nseeds),
    warping.method = list(warping_method),
    similarity.method = list(similarity_method),
    center.method = list(center_method),
    out
  )

  class(out) <- "kmap"

  out
}

#'@title Clustering and alignment of functional data
#'
#'@description kmap jointly performs clustering and alignment of a functional dataset
#'(multidimensional or unidimensional functions). To run kmap function with different
#'numbers of clusters and/or different alignment methods see
#'
#'@usage res<-kmap(  x, y,  y1, n_clust,warping_method, center_method,similarity_method,
#'optim_method, seeds, span, delta, d_max, s_max, n_out, toll, fence,
#' iter_max,show_iter, check_total_similarity)
#'
#'@param x A matrix \emph{n.func} X \emph{grid.size} or vector \emph{grid.size}:
#'     the abscissa values where each function is evaluated. \emph{n.func}: number
#'     of functions in the dataset. \emph{grid.size}: maximal number of abscissa values
#'      where each function is evaluated. The abscissa points may be unevenly spaced and
#'      they may differ from function to function. \code{x} can also be a vector of length
#'       \emph{grid.size}. In this case, \code{x} will be used as abscissa grid for all
#'       functions.Furthermore if the grid's size differs from one function to another the
#'       matrix must be completed with NA values.The parameter \code{y} must be
#'       provided.
#'
#'@param y A matrix \emph{n.func} X \emph{grid.size} or array \emph{n.func} X \emph{grid.size}
#'         X \emph{d}: evaluations of the set of original functions on the abscissa grid
#'         \code{x}. \emph{n.func}: number of functions in the dataset. \emph{grid.size}:
#'         maximal number of abscissa values where each function is evaluated. \emph{d}:
#'         (only if the sample is multidimensional) number of function components, i.e.
#'          each function is a \emph{d}-dimensional curve. The parameter \code{yx} must be
#'           provided.
#'
#' @param seeds vector max(n.clust) or matrix nstart X n.clust: indexes of the functions
#'              to be used as initial centers. If it is a matrix, each row contains the
#'               indexes of the initial centers of one of the nstart initializations.
#'               In the case where not all the values of seeds are provided, those not
#'               provided are randomly chosen among the n.func original functions.
#'               If seeds=NULL all the centers are randomly chosen. Default value of seeds
#'               is NULL.
#'
#'@param n_clust scalar: required number of clusters. Default value is 1. Note that if
#'               n.clust=1 kma performs only alignment without clustering.
#'
#' @param warping_method character: type of alignment required. If warping.method='noalign'
#'                       kma performs only clustering (without alignment).
#'                       If warping.method='affine' kma performs alignment
#'                       (and possibly clustering) of functions using linear affine
#'                       transformation as warping functions, i.e., x.final = dilation*x + shift.
#'                        If warping.method='shift' kma allows only shift, i.e.,
#'                        x.final = x + shift. If warping.method='dilation' kma allows only
#'                        dilation, i.e., x.final = dilation*x. Default value is 'affine'.
#'@param center_method character: type of clustering method to be used.
#'                     Possible choices are: 'mean' and 'medoid' and 'pseudomedoid'.
#'                     Default value is 'mean'.
#'
#'@param similarity_method character: required similarity measure. Possible choices are:
#'                         'pearson','l2'. Default value is 'pearson'.
#'                         See (what to define?) for details.
#'@param optim_method character: optimization method chosen to find the best warping
#'                    functions at each iteration. Possible choices are: 'bobyqa'.
#'                     See optim function for details.
#'                     Default method is 'bobyqa'
#'
#'@param warping_opt numeric vector. The parameters to set depend on the warping_method chosen.
#'                   If warping_method ='affine' warping_opt <- c( max_dilation , max_shift).
#'                   If warping_method <- 'dilation' warping_opt <- c(max_dilation).
#'                   If warping_method <- 'shift'  warping_opt <- c(max_shift).
#'                   If warping_method <- 'noalign' warping_opt <- as.numeric().
#'
#'@param center_opt numeric vector. The parameters to set depend on the center_method chosen.
#'                  If center_method ='mean' center_opt <- c( span, delta).
#'                  If center_method ='medoid' center_opt <- as.numeric().
#'                  If center_method ='pseudomedoid' center_opt <- as.numeric().
#'
#'@param out_opt numeric vector. The parameters to set are (n_out , tollerance, max_iteration).
#'               n_out is the size of the grid where the centers will be computed.
#'               tollerance is a stop condition parameter.
#'               max_iterationa is a stop condition parameter.
#'
#'@param fence boolean: if fence=TRUE a control is activated at the end of each iteration.
#'             The aim of the control is to avoid warping outliers with respect
#'              to their computed distributions. If fence=TRUE the running time can increase
#'               considerably. Default value of fence is FALSE.
#'
#'@param check_total_similarity boolean: if check.total.similarity=TRUE at each iteration
#'                              the algorithm checks if there is a decrease of the total
#'                              similarity and stops. In the affermative case the result
#'                              obtained in the penultimate iteration is returned.
#'                               Defaultvalue is TRUE.
#'
#'@param show_iter boolean: if show.iter=TRUE kmap shows the current iteration of the
#'                 algorithm. Default value is TRUE.
#'
#'@param comp_original_center boolean: if comp_original_center=TRUE the initial center with relative dissimilarities
#'                             is computed otherwise this step is skipped. Default value is FALSE.
#'@param par_opt  numeric vector. Parallel options. The parameters to set are (num_threads,parallel_version)
#'                parallel_version  available are 0 and 1 :  0 is a trivial parallelization in which each
#'                thread compute the center of a cluster; 1 is a more efficient parallelization in which
#'                all the threads compute the centers sequentially(available only with center_method = 'medoid').
#'
#'@return The function output is a list containing the following elements:
#'
#' \item{ x }{ as input.}
#' \item{ y }{ as input. }
#' \item{seeds}{ vector with the indeces used in the algorithm.}
#'
#' \item{ warping.method }{ as input.}
#' \item{ similarity.method }{ as input. }
#' \item{ center.method }{ as input. }
#' \item{ iterations }{scalar: total number of iterations performed by kma function.}
#' \item{ n.clust }{ as input. }
#' \item{ x.center.orig }{numeric vector \emph{n_out}: abscissa of the center computed if \emph{comp_original_center}=TRUE.}
#' \item{ y.center.orig }{numeric vector \emph{n_out} or matrix \emph{n_out} X \emph{n_dim}: value of the center computed if \emph{comp_original_center}=TRUE.}
#' \item{similarity.origin}{numeric vector \emph{n_obs} dissimilarity,similarity or distance of the original center respect the obserbations computed if \emph{comp_original_center}=TRUE.}
#' \item{ x.final }{matrix n.func X grid.size: aligned abscissas.}
#' \item{ n.clust.final }{ scalar: final number of clusters. Note that n.clust.final may differ from initial number of clusters (i.e.,from n.clust) if some clusters are found to be empty.}
#' \item{ x.centers.final }{matrix n.clust.final X grid.size: abscissas of the final function centers.}
#' \item{ y.centers.final }{matrix n.clust.final X n.out or array n.clust.final X n.out x n_dim , contain the evaluations of the final functions centers.}
#' \item{ templates_vec }{list iteration : each element of the list contain centers of that iteration.}
#' \item{ x_out_vec }{list iteration : each element of the list contain the abscissa of the centers of that iteration.}
#' \item{ labels }{vector n_obs: cluster assignments.}
#' \item{ similarity.final }{vector n_obs: similarities,dissimilarities or distance between each function and the center of the cluster the function is assigned to.}
#' \item{parameters.list}{list iterations: warping parameters at each iteration.}
#' \item{parameters}{matrix n_par X n_obs: warping parameters applied to the original abscissas x to obtain the aligned abscissas x.final.}
#' \item{timer}{vector : time of execution. }


kmap <- function(x, y, seeds= NULL, n_clust = 1,
       warping_method ='affine', center_method ='mean',
       similarity_method ='pearson', optim_method = 'bobyqa',
       warping_opt=c(0.15,0.15), center_opt = c(0.01,0.1), out_opt = c(100 , 0.001 , 100),
       fence = FALSE, check_total_similarity = TRUE,show_iter = TRUE,
       comp_original_center=FALSE, par_opt=c(1,0))
{

  if(is.null(y))
    stop("Provvide a valid function")

  if(length(dim(y))==2){
    y<-array(y,c(dim(y)[1],dim(y)[2],1))
  }

  if(is.vector(x))
   x <- t(replicate(dim(y)[1],x,1))

  if(is.null(seeds)){
    nseeds<-sample(0:(nrow(y)-1),n_clust)
  }else{
    nseeds<-seeds
  }
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

  out<-c(x=list(x),y=list(y), seeds=list(nseeds), warping.method = list(warping_method),
         similarity.method = list(similarity_method),center.method = list(center_method),out)

}

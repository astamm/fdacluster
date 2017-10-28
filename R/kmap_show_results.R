#'@title plot results kmap
#'
#'@description Show results of the clustering with alignment of functional data
#
#'@usage kma_show_results <-function (Result, bp_sim)
#'@param Result output of kmap.
#'
#'@param bp_sim boolean: if TRUE dissimilarity,similarity or distance boxplot are plotted.
#'

kmap_show_results <-function (Result, bp_sim=FALSE){
  # Result<-res
  lwd.functions = 1
  lwd.centers = 3


# settings ----------------------------------------------------------------
  if (length(Result) == 0) {
    stop("First parameter (Result) must not be NULL")
  }

    iterations <- Result$iterations
    x <- as.matrix(Result$x)
    y0 <- Result$y
    n.clust <- Result$n.clust

    x.center.orig <- Result$x.center.orig
    y0.center.orig <- Result$y.center.orig
    sims.orig <- Result$similarity.orig
    cen_or<- TRUE

    if(length(x.center.orig)==0)
    cen_or<- FALSE

    x.final <- Result$x.final
    n.clust.final <- Result$n.clust.final
    x.centers.final <- Result$x.centers.final
    y0.centers.final <- Result$y.centers.final
    labels <- Result$labels
    sims.final <- Result$similarity.final
    parameters.list <- Result$parameters.list
    parameters <- Result$parameters
    warping.method <- Result$warping.method
    similarity.method <- Result$similarity.method
    center.method <- Result$center.method



    #
    # LEGGO DIMENSIONI
    #

    sim.final <- mean(sims.final)
    labels.unique <- sort(unique(labels))

    n.obs <- dim(y0)[1]
    n.camp <- dim(y0)[2]
    n.dim <- dim(y0)[3]



    ## COLORI
    myrainbow <- c("red", "blue", "green3", "orange", "grey", "yellow")
    myrainbow.dark <- c("darkred", "darkblue", "darkgreen", "darkorange","black", "brown")
    colori.dopo <- rainbow(length(labels.unique))
    myrainbow <- c(myrainbow, colori.dopo)
    myrainbow <- myrainbow[1:n.clust.final]
    myrainbow.dark <- c(myrainbow.dark, colori.dopo)
    myrainbow.dark <- myrainbow.dark[1:n.clust.final]
    colours.random <- rainbow(n.obs)

    colours.bygroup <- rep(0, n.obs)
    colours.bygroup.dark <- rep(0, n.obs)

    for (k in labels.unique) {
        colours.bygroup[which(labels == k)] <- myrainbow[k]
        colours.bygroup.dark[which(labels == k)] <- myrainbow.dark[k]
    }

    colours.templates.iter1 <- myrainbow
    colours.templates.last <- myrainbow
    colours.warping <- rep(0, n.obs)

    for (k in labels.unique) {
        colours.warping[which(labels == k)] <- myrainbow[k]
        colours.templates.iter1[k] <- myrainbow[k]
        colours.templates.last[k] <- colours.bygroup[which(labels == k)[1]]
    }




# print original data -----------------------------------------------------

    dev.new()
    par(mfrow=c(2,n.dim))

    for (l in 1:n.dim) {
      var<-c("x","y","z")

      matplot(t(x), t(y0[, , l]), type = "l", lwd = lwd.functions,
              col = colours.random,xlab = "t", ylab = var[l])

      title(main = paste(c("Original Data, dimension ", l), collapse = ""))
      if(cen_or)
        lines(x.center.orig, y0.center.orig[l, ],lwd = lwd.centers, col = colours.templates.iter1)
    }



# printed aligned data ----------------------------------------------------

    for (l in 1:n.dim) {

      tex <- paste(c("k = ", n.clust), collapse = "")
      matplot(t(x.final), t(y0[, , l]), type = "l",lwd = lwd.functions, col = colours.bygroup.dark,
                  xlab = "x", ylab = tex)

        title2 <- c("Registration: ", warping.method)
      title2 <- paste(title2, collapse = "")

      title22 <- paste(c("Aligned Data; dimension ", l), collapse = "")
      title2def <- c(title2, title22)
      title(main = title2def)

      for (k in labels.unique) {
        col_cl<-myrainbow
        lines(x.centers.final, y0.centers.final[k,,l], lwd = lwd.centers, col = col_cl[k])
      }

      text <- rep(0, length(labels.unique))
      for (i in 1:length(labels.unique)) {
        tt <- c("Cluster ", i)
        tt <- paste(tt, collapse = "")
        text[i] <- tt
      }
      lty <- rep(1, length(text))
      legend("topleft", legend = text, col = colours.templates.last,
             lty = lty, cex = 0.6)
    }








# print warping functions -------------------------------------------------


# #
# #      print warping functions
# #
#   dev.new()
#   plot(t(x[1, ]), t(x.final[1, ]), xlim = c(min(x, na.rm = TRUE),
#       max(x, na.rm = TRUE)), ylim = c(min(x, na.rm = TRUE),
#       max(x, na.rm = TRUE)), type = "l", lwd = 1, col = colours.warping[1],
#       xlab = "x", ylab = "y", asp = 1)
#     for (i in 2:n.obs) {
#         lines(t(x[i, ]), t(x.final[i, ]), type = "l", lwd = 1,
#             col = colours.warping[i], xlab = "x", ylab = "y",
#             asp = 1)
#     }
#     title3 <- c("Registration: ", Result$warping.method)
#     title3 <- paste(title3, collapse = "")
#     title33 <- c("Warping Functions")
#     title3def <- c(title3, title33)
#     title(main = title3def)
#     abline(v = min(x))
#     abline(v = max(x))


# print similarity boxes --------------------------------------------------
    #
    #  print similarty boxes
    #
   if(bp_sim == TRUE){
    if(cen_or==TRUE){

      dev.new()

      if (similarity.method == "pearson") {
        boxplot(-sims.orig, -sims.final, notch = FALSE, boxwex = 0.3,
                col = c("grey", "orange"), ylim = c(min(-sims.orig,
                                                        -sims.final), 1))
      }
      if (similarity.method == "l2") {
        boxplot(sims.orig, sims.final, notch = FALSE, boxwex = 0.3,
                col = c("grey", "orange"), ylim = c(min(sims.orig,sims.final),
                                                    max(sims.orig, sims.final)))
      }


     # creo testo
      text <- rep(0, length(labels.unique))
      for (i in 1:length(labels.unique)) {
        tt <- c("Cluster ", i)
        tt <- paste(tt, collapse = "")
        text[i] <- tt
      }

      lty <- rep(1, length(text))
      legend(x = min(x), y = max(x), legend = text, col = colours.templates.last,
      lty = lty, cex = 0.6)
      etichette <- rep(0, 2)
      etichette[1] <- "orig. data"
      etichette[2] <- paste("k =", n.clust)
      title4 <- c("Registration: ", Result$warping.method)
      title4 <- paste(title4, collapse = "")
      title44 <- c("Boxplot Similarity Indexes")
      title4def <- c(title4, title44)
      title(main = title4def)
      axis(1, at = 1:2, labels = etichette, las = 0)
    }else{
      print("Sim_boxes: Original center not available")
    }
   }


}

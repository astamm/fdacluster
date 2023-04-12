library(fdacluster)
library(future)
library(progressr)
handlers("rstudio")
ncores <- max(parallel::detectCores() - 1L, 1L)
plan(multisession, workers = ncores)
with_progress({
  sim30_mcaps <- compare_caps(
    x = simulated30$x,
    y = simulated30$y,
    warping_class = c("none", "shift", "dilation", "affine"),
    clustering_method = "kmeans",
    centroid_type = "mean"
  )
})
plan(sequential)
usethis::use_data(sim30_mcaps, overwrite = TRUE, compress = "xz", version = 3)

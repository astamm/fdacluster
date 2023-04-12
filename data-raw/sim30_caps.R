library(fdacluster)
library(future)
library(progressr)
handlers("rstudio")
ncores <- max(parallel::detectCores() - 1L, 1L)
plan(multisession, workers = ncores)
with_progress({
  sim30_caps <- fdakmeans(
    x = simulated30$x,
    y = simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson"
  )
})
plan(sequential)
usethis::use_data(sim30_caps, overwrite = TRUE, compress = "xz", version = 3)

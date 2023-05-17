library(fdacluster)
library(future)
ncores <- max(parallel::detectCores() - 1L, 1L)

plan(multisession, workers = ncores)

amplitude_data <- compare_caps(
  x = simulated30_sub$x,
  y = simulated30_sub$y,
  n_clusters = 1:5,
  metric = "l2",
  warping_class = c("none", "shift", "dilation", "affine"),
  clustering_method = "hclust-complete",
  centroid_type = "mean",
  cluster_on_phase = FALSE
)

phase_data <- compare_caps(
  x = simulated30_sub$x,
  y = simulated30_sub$y,
  n_clusters = 1:5,
  metric = "l2",
  warping_class = c("shift", "dilation", "affine"),
  clustering_method = "hclust-complete",
  centroid_type = "mean",
  cluster_on_phase = TRUE
)

plan(sequential)

usethis::use_data(
  amplitude_data,
  phase_data,
  internal = TRUE,
  overwrite = TRUE,
  compress = "xz",
  version = 3
)

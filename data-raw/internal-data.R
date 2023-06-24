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

out_manual <- fdakmeans(
  x = simulated30$x,
  y = simulated30$y,
  n_clusters = 2,
  seeds = c(1, 21),
  warping_class = "affine",
  centroid_type = "mean",
  metric = "l2",
  cluster_on_phase = FALSE,
  use_verbose = FALSE
)

withr::with_seed(1234, {
  initial_seeds <- replicate(10, sample.int(30, 2, replace = FALSE), simplify = FALSE)
  outs_manual <- purrr::map(initial_seeds, \(.seeds) {
    fdakmeans(
      x = simulated30$x,
      y = simulated30$y,
      n_clusters = 2,
      seeds = .seeds,
      warping_class = "affine",
      centroid_type = "mean",
      metric = "l2",
      cluster_on_phase = FALSE,
      use_verbose = FALSE
    )
  })
})

withr::with_seed(1234, {
  outs_kpp <- replicate(10, {
    fdakmeans(
      x = simulated30$x,
      y = simulated30$y,
      n_clusters = 2,
      seeding_strategy = "kmeans++",
      warping_class = "affine",
      centroid_type = "mean",
      metric = "l2",
      cluster_on_phase = FALSE,
      use_verbose = FALSE
    )
  }, simplify = FALSE)
})

withr::with_seed(1234, {
  outs_ekpp <- replicate(10, {
    fdakmeans(
      x = simulated30$x,
      y = simulated30$y,
      n_clusters = 2,
      seeding_strategy = "exhaustive-kmeans++",
      warping_class = "affine",
      centroid_type = "mean",
      metric = "l2",
      cluster_on_phase = FALSE,
      use_verbose = FALSE
    )
  }, simplify = FALSE)
})

out_hclust <- fdakmeans(
  x = simulated30$x,
  y = simulated30$y,
  n_clusters = 2,
  seeding_strategy = "hclust",
  warping_class = "affine",
  centroid_type = "mean",
  metric = "l2",
  cluster_on_phase = FALSE,
  use_verbose = FALSE
)

plan(sequential)

usethis::use_data(
  amplitude_data,
  phase_data,
  initial_seeds,
  out_manual,
  outs_manual,
  outs_kpp,
  outs_ekpp,
  out_hclust,
  internal = TRUE,
  overwrite = TRUE,
  compress = "xz",
  version = 3
)

library(fdacluster)
library(future)
ncores <- max(parallel::detectCores() - 1L, 1L)

plan(multisession, workers = ncores)

# hclust ------------------------------------------------------------------

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

# kmeans-initialisation ---------------------------------------------------

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

# berkeley-growth ---------------------------------------------------------

growth <- fda::growth
mb <- as.factor(c(rep("male", dim(growth$hgtm)[2]), rep("female", dim(growth$hgtf)[2])))
N <- length(mb)
x <- growth$age
M <- length(x)
y0 <- cbind(growth$hgtm, growth$hgtf)
basisobj <- fda::create.bspline.basis(rangeval = range(x), nbasis = 15)
fd_vals <- purrr::map(1:N, \(n) {
  yobs <- y0[, n]
  result <- fda::smooth.basis(x, yobs, basisobj)
  yfd <- result$fd
  cost <- function(lam) {
    yfdPar <- fda::fdPar(yfd, 2, lam)
    out <- fda::smooth.basis(x, yobs, yfdPar)
    out$gcv
  }
  lambda_opt <- stats::optimise(cost, c(1e-8, 1))$minimum
  if (lambda_opt <= 1e-8)
    cli::cli_alert_warning("The optimal penalty has reached the lower bound (1e-8) for curve #{n}.")
  if (lambda_opt >= 1)
    cli::cli_alert_warning("The optimal penalty has reached the upper bound (1) for curve #{n}.")
  yfdPar <- fda::fdPar(yfd, 2, lambda_opt)
  fda::smooth.fd(yfd, yfdPar)
})
fd <- fda::fd(
  coef = fd_vals |>
    purrr::map("coefs") |>
    purrr::reduce(cbind),
  basisobj = basisobj
)
y0 <- fda::eval.fd(x, fd, 0)
y1 <- fda::eval.fd(x, fd, 1)

growth_mcaps <- compare_caps(
  x = x,
  y = t(y1),
  n_clusters = 2,
  metric = "l2",
  centroid_type = "mean",
  warping_class = "affine",
  cluster_on_phase = TRUE
)

growth_caps <- fdahclust(
  x = x,
  y = t(y1),
  n_clusters = 2,
  metric = "l2",
  warping_class = "affine",
  centroid_type = "mean",
  cluster_on_phase = TRUE
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
  growth_mcaps,
  growth_caps,
  internal = TRUE,
  overwrite = TRUE,
  compress = "xz",
  version = 3
)

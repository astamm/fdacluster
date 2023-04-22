test_that('`fdakmeans()` works with fda::fd input object.', {
  fd <- fda::as.fd(fda::smooth.basisPar(
    simulated30_sub$x[1, ],
    t(simulated30_sub$y[, 1, ]),
    lambda = 0.00001)
  )
  out <- fdakmeans(
    x = simulated30_sub$x,
    y = fd,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with funData::funData input object.', {
  fd <- funData::funData(simulated30_sub$x[1, ], simulated30_sub$y[, 1, ])
  out <- fdakmeans(
    x = fd,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with fixed initial seeds.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with kmeans++ seeding strategy.', {
  withr::with_seed(1234, {
    out <- fdakmeans(
      simulated30_sub$x,
      simulated30_sub$y,
      seeding_strategy = "kmeans++",
      n_clusters = 2,
      centroid_type = "mean",
      warping_class = "affine",
      metric = "pearson",
      use_verbose = FALSE
    )
  })

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with exhaustive-kmeans++ seeding strategy.', {
  library(future)
  ncores <- max(parallel::detectCores() - 1, 1)
  plan(multisession, workers = ncores)
  withr::with_seed(1234, {
    out <- fdakmeans(
      simulated30_sub$x,
      simulated30_sub$y,
      seeding_strategy = "exhaustive-kmeans++",
      n_clusters = 2,
      centroid_type = "mean",
      warping_class = "affine",
      metric = "pearson",
      use_verbose = FALSE
    )
  })
  plan(sequential)

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with exhaustive seeding strategy.', {
  library(future)
  ncores <- max(parallel::detectCores() - 1, 1)
  plan(multisession, workers = ncores)
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeding_strategy = "exhaustive",
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )
  plan(sequential)

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 3)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with hclust seeding strategy.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeding_strategy = "hclust",
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with dilation warping.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "dilation",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with no warping.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "none",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 4)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with shift warping.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "shift",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with srsf warping.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "srsf",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 1)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with medoid centroid.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with lowess centroid.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "lowess",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with poly centroid.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "poly",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 1)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with l2 metric.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "l2",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works when clustering on phase.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 11, 21),
    n_clusters = 3,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    cluster_on_phase = TRUE,
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(3, 1, 50))
  expect_equal(dim(out$grids), c(3, 50))
  expect_equal(out$n_clusters, 3)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works in verbose mode.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = TRUE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with parallel mode on distance calculation.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE,
    parallel_method = 1L
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with fence adaptive algorithm.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE,
    use_fence = TRUE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with computation of overall center.', {
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE,
    compute_overall_center = TRUE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "warpings", "grids", "n_clusters", "memberships",
                      "distances_to_center", "silhouettes",
                      "amplitude_variation", "total_variation", "n_iterations",
                      "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 50))
  expect_equal(dim(out$aligned_curves), c(30, 1, 50))
  expect_equal(dim(out$center_curves), c(2, 1, 50))
  expect_equal(dim(out$grids), c(2, 50))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 50))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

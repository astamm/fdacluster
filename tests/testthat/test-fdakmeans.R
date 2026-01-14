test_that('`fdakmeans()` works with fda::fd input object.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  fd <- fda::as.fd(fda::smooth.basisPar(
    simulated30_sub$x[1, ],
    t(simulated30_sub$y[, 1, ]),
    lambda = 0.00001
  ))
  out <- fdakmeans(
    x = simulated30_sub$x,
    y = fd,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with funData::funData input object.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  fd <- funData::funData(simulated30_sub$x[1, ], simulated30_sub$y[, 1, ])
  out <- fdakmeans(
    x = fd,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with fixed initial seeds.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with kmeans++ seeding strategy.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  withr::with_seed(1234, {
    out <- fdakmeans(
      simulated30_sub$x,
      simulated30_sub$y,
      seeding_strategy = "kmeans++",
      n_clusters = K,
      centroid_type = "mean",
      warping_class = "affine",
      metric = "pearson",
      use_verbose = FALSE
    )
  })

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with exhaustive-kmeans++ seeding strategy.', {
  skip_if_not_installed("lpSolve")
  skip_on_cran()

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  library(future)
  ncores <- max(parallel::detectCores() - 1, 1)
  plan(multisession, workers = ncores)
  withr::with_seed(1234, {
    out <- fdakmeans(
      simulated30_sub$x,
      simulated30_sub$y,
      seeding_strategy = "exhaustive-kmeans++",
      n_clusters = K,
      centroid_type = "mean",
      warping_class = "affine",
      metric = "pearson",
      use_verbose = FALSE
    )
  })
  plan(sequential)

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with exhaustive seeding strategy.', {
  skip_if_not_installed("lpSolve")
  skip_on_cran()

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  library(future)
  ncores <- max(parallel::detectCores() - 1, 1)
  plan(multisession, workers = ncores)
  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeding_strategy = "exhaustive",
    n_clusters = K,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )
  plan(sequential)

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with hclust seeding strategy.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeding_strategy = "hclust",
    n_clusters = K,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with dilation warping.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "mean",
    warping_class = "dilation",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with no warping.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "mean",
    warping_class = "none",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 4)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with shift warping.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "mean",
    warping_class = "shift",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with boundary-preserving diffeomorphism warping.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = K,
    is_domain_interval = TRUE,
    transformation = "srvf",
    centroid_type = "mean",
    warping_class = "bpd",
    metric = "l2",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 1)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with median centroid.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "median",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with medoid centroid.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "medoid",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with lowess centroid.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "lowess",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with poly centroid.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "poly",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with normalized l2 metric.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "normalized_l2",
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works when clustering on phase.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 3L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 11, 21),
    n_clusters = K,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    cluster_on_phase = TRUE,
    use_verbose = FALSE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works in verbose mode.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = TRUE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with parallel mode on distance calculation.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "medoid",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE,
    parallel_method = 1L
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with fence adaptive algorithm.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE,
    use_fence = TRUE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with computation of overall center.', {
  skip_if_not_installed("lpSolve")

  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = K,
    centroid_type = "mean",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE,
    compute_overall_center = TRUE
  )

  expect_true(is_caps(out))
  expect_equal(length(out), 14)
  expected_names <- c(
    "original_curves",
    "original_grids",
    "aligned_grids",
    "center_curves",
    "center_grids",
    "n_clusters",
    "memberships",
    "distances_to_center",
    "silhouettes",
    "amplitude_variation",
    "total_variation",
    "n_iterations",
    "call_name",
    "call_args"
  )
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$original_grids), c(N, P))
  expect_equal(dim(out$aligned_grids), c(N, P))
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$center_grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that("`fdahclust()` works", {
  dims <- dim(simulated30_sub$y)
  N <- dims[1]
  L <- dims[2]
  P <- dims[3]
  K <- 2L
  out <- fdahclust(
    x = simulated30_sub$x,
    y = simulated30_sub$y,
    n_clusters = K,
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
  expect_equal(dim(out$original_curves), dims)
  expect_equal(dim(out$aligned_curves), dims)
  expect_equal(dim(out$center_curves), c(K, L, P))
  expect_equal(dim(out$grids), c(K, P))
  expect_equal(out$n_clusters, K)
  expect_equal(length(out$memberships), N)
  expect_equal(length(out$distances_to_center), N)
  expect_equal(dim(out$warpings), c(N, P))
  expect_equal(out$n_iterations, 0)
  expect_equal(out$call_name, "fdahclust")
  expect_true(inherits(out$call_args, "list"))
})

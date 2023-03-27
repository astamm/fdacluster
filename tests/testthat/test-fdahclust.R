test_that("`fdahclust()` works", {
  out <- fdahclust(
    x = simulated30$x,
    y = simulated30$y,
    n_clusters = 2,
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
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(dim(out$grids), c(2, 200))
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 0)
  expect_equal(out$call_name, "fdahclust")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "medoid"`,
          `warping_class = "affine"` and `distance = "pearson"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "affine",
    distance = "pearson"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "medoid"`,
          `warping_class = "affine"` and `distance = "l2"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "affine",
    distance = "l2"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "medoid"`,
          `warping_class = "dilation"` and `distance = "pearson"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "dilation",
    distance = "pearson"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "medoid"`,
          `warping_class = "dilation"` and `distance = "l2"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "dilation",
    distance = "l2"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "medoid"`,
          `warping_class = "shift"` and `distance = "pearson"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "shift",
    distance = "pearson"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "medoid"`,
          `warping_class = "shift"` and `distance = "l2"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "shift",
    distance = "l2"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "medoid"`,
          `warping_class = "none"` and `distance = "pearson"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "affine",
    distance = "pearson"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "medoid"`,
          `warping_class = "none"` and `distance = "l2"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "affine",
    distance = "pearson"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "mean"`,
          `warping_class = "affine"` and `distance = "pearson"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "mean",
    warping_class = "affine",
    distance = "pearson"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "mean"`,
          `warping_class = "affine"` and `distance = "l2"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "affine",
    distance = "l2"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "mean"`,
          `warping_class = "dilation"` and `distance = "pearson"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "dilation",
    distance = "pearson"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "mean"`,
          `warping_class = "dilation"` and `distance = "l2"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "dilation",
    distance = "l2"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "mean"`,
          `warping_class = "shift"` and `distance = "pearson"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "shift",
    distance = "pearson"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "mean"`,
          `warping_class = "shift"` and `distance = "l2"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "shift",
    distance = "l2"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "mean"`,
          `warping_class = "none"` and `distance = "pearson"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "affine",
    distance = "pearson"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `centroid_type = "mean"`,
          `warping_class = "none"` and `distance = "l2"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "affine",
    distance = "pearson"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 2)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

test_that('`fdakmeans()` works with `warping_class = "srsf"`.', {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    warping_class = "srsf"
  )

  expect_true(is_kmaps(out))
  expect_equal(length(out), 11)
  expected_names <- c("original_curves", "aligned_curves", "center_curves",
                      "grid", "n_clusters", "memberships", "distances_to_center",
                      "warpings", "n_iterations", "call_name", "call_args")
  expect_equal(names(out), expected_names)
  expect_equal(dim(out$original_curves), c(30, 1, 200))
  expect_equal(dim(out$aligned_curves), c(30, 1, 200))
  expect_equal(dim(out$center_curves), c(2, 1, 200))
  expect_equal(length(out$grid), 200)
  expect_equal(out$n_clusters, 2)
  expect_equal(length(out$memberships), 30)
  expect_equal(length(out$distances_to_center), 30)
  expect_equal(dim(out$warpings), c(30, 200))
  expect_equal(out$n_iterations, 3)
  expect_equal(out$call_name, "fdakmeans")
  expect_true(inherits(out$call_args, "list"))
})

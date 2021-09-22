library(fdacluster)
context("Warping Method")

test_that(" the warping methods work", {
  expect_equal(length(
    kma(
      x = aneurisk65$x,
      y = aneurisk65$y,
      seeds = NULL,
      n_clust = 2,
      warping_method = "affine",
      use_verbose = FALSE
    )
  ), 23)

  expect_equal(length(
    kma(
      x = aneurisk65$x,
      y = aneurisk65$y,
      seeds = NULL,
      n_clust = 2,
      warping_method = "dilation",
      use_verbose = FALSE
    )
  ), 23)

  expect_equal(length(
    kma(
      x = aneurisk65$x,
      y = aneurisk65$y,
      seeds = NULL,
      n_clust = 2,
      warping_method = "shift",
      use_verbose = FALSE
    )
  ), 23)

  expect_equal(length(
    kma(
      x = aneurisk65$x,
      y = aneurisk65$y,
      seeds = NULL,
      n_clust = 2,
      warping_method = "none",
      use_verbose = FALSE
    )
  ), 23)
})

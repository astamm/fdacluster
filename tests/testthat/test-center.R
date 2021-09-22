library(fdacluster)
context("Center Method")

test_that(" the warping methods work", {
  expect_equal(length(
    kma(
      x = aneurisk65$x,
      y = aneurisk65$y,
      seeds = NULL,
      n_clust = 2,
      center_method = "mean",
      use_verbose = FALSE
    )
  ), 23)

  expect_equal(length(
    kma(
      x = aneurisk65$x,
      y = aneurisk65$y,
      seeds = NULL,
      n_clust = 2,
      center_method = "medoid",
      use_verbose = FALSE
    )
  ), 23)
})

library(fdakmapp)
context("Warping Method")

test_that(" the warping methods work", {
  expect_equal(length(
    kmap(
      x = aneurisk65$x,
      y = aneurisk65$y,
      seeds = NULL,
      n_clust = 2,
      warping_method = "affine",
      warping_opt = c(0.15, 0.15),
      show_iter = FALSE
    )
  ), 22)

  expect_equal(length(
    kmap(
      x = aneurisk65$x,
      y = aneurisk65$y,
      seeds = NULL,
      n_clust = 2,
      warping_method = "dilation",
      warping_opt = c(0.15),
      show_iter = FALSE
    )
  ), 22)

  expect_equal(length(
    kmap(
      x = aneurisk65$x,
      y = aneurisk65$y,
      seeds = NULL,
      n_clust = 2,
      warping_method = "shift",
      warping_opt = c(0.15),
      show_iter = FALSE
    )
  ), 22)

  expect_equal(length(
    kmap(
      x = aneurisk65$x,
      y = aneurisk65$y,
      seeds = NULL,
      n_clust = 2,
      warping_method = "noalign",
      warping_opt = as.numeric(),
      show_iter = FALSE
    )
  ), 22)
})

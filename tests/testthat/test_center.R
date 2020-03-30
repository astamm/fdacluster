library(fdakmapp)
context("Center Method")

test_that(" the warping methods work", {
  expect_equal(length(
    kmap(
      x = aneurisk65$x,
      y = aneurisk65$y,
      seeds = NULL,
      n_clust = 2,
      center_method = "mean",
      center_opt = c(0.01, 0.1),
      show_iter = FALSE
    )
  ), 22)

  expect_equal(length(
    kmap(
      x = aneurisk65$x,
      y = aneurisk65$y,
      seeds = NULL,
      n_clust = 2,
      center_method = "medoid",
      center_opt = as.numeric(),
      show_iter = FALSE
    )
  ), 22)

  expect_equal(length(
    kmap(
      x = aneurisk65$x,
      y = aneurisk65$y,
      seeds = NULL,
      n_clust = 2,
      center_method = "medoid",
      center_opt = as.numeric(),
      show_iter = FALSE,
      par_opt = c(4, 1)
    )
  ), 22)
})

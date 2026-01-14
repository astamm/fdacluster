test_that("`fdadist()` works", {
  skip_if_not_installed("lpSolve")

  D <- fdadist(simulated30_sub$x, simulated30_sub$y)
  expect_true(inherits(D, "dist"))
  expect_equal(length(D), 435)
})

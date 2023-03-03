test_that("`fdadist()` works", {
  D <- fdadist(simulated30$x, simulated30$y)
  expect_true(inherits(D, "dist"))
  expect_equal(length(D), 435)
})

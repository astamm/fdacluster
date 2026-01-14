test_that("`diagnostic_plot()` works", {
  skip_if_not_installed("lpSolve")

  out <- fdakmeans(
    simulated30_sub$x,
    simulated30_sub$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "affine",
    metric = "pearson",
    use_verbose = FALSE
  )
  p <- diagnostic_plot(out)
  expect_equal(dim(p$data), c(60, 4))
  expect_equal(names(p$data), c("Membership", "ObservationID", "name", "value"))
})

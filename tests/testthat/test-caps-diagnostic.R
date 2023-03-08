test_that("`diagnostic_plot()` works", {
  out <- fdakmeans(
    simulated30$x,
    simulated30$y,
    seeds = c(1, 21),
    n_clusters = 2,
    centroid_type = "medoid",
    warping_class = "affine",
    metric = "pearson"
  )
  p <- diagnostic_plot(out)
  expect_equal(dim(p$data), c(60, 4))
  expect_equal(names(p$data), c("Membership", "ObservationID", "name", "value"))
})

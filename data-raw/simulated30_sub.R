N <- nrow(simulated30$x)
P <- 50L
x_in <- simulated30$x
x_out <- matrix(nrow = N, ncol = P)
y_in <- simulated30$y
y_out <- array(dim = c(N, 1L, P))
for (n in 1:N) {
  x_out[n, ] <- seq(min(x_in[n, ]), max(x_in[n, ]), length.out = P)
  y_out[n, 1, ] <- approx(x_in[n, ], y_in[n, 1, ], xout = x_out[n, ])$y
}
simulated30_sub <- NULL
simulated30_sub$x <- x_out
simulated30_sub$y <- y_out
usethis::use_data(simulated30_sub, overwrite = TRUE, compress = "xz", version = 3)

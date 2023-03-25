impute_via_mean <- function(X, labels) {
  N <- nrow(X)
  M <- ncol(X)
  idx_min <- rep(NA, N)
  idx_max <- rep(NA, N)

  for (n in 1:N) {
    break_low <- FALSE
    break_high <- FALSE
    for (m in 1:M) {
      if (!is.na(X[n, m]) && !break_low) {
        idx_min[n] <- m
        break_low <- TRUE
      }
      if (!is.na(X[n, M - m + 1]) && !break_high) {
        idx_max[n] <- M - m + 1
        break_high <- TRUE
      }
      if (break_low && break_high)
        break
    }
  }

  Xout <- X
  prior <- matrix(nrow = 0, ncol = 4)

  for (m in (max(idx_min) + 1):M) {
    idx_na <- which(is.na(X[, m]))
    for (n in idx_na) {
      subs <- which(labels == labels[n])
      spl <- X[subs, m] - X[subs, m - 1]
      Xout[n, m] <- mean(spl, na.rm = TRUE) + Xout[n, m - 1]
      prior <- rbind(prior, c(n, m, Xout[n, m], stats::sd(spl, na.rm = TRUE)))
    }
  }

  for (m in (min(idx_max) - 1):1) {
    idx_na <- which(is.na(X[, m]))
    for (n in idx_na) {
      subs <- which(labels == labels[n])
      spl <- X[subs, m] - X[subs, m + 1]
      Xout[n, m] <- mean(spl, na.rm = TRUE) + Xout[n, m + 1]
      prior <- rbind(prior, c(n, m, Xout[n, m], stats::sd(spl, na.rm = TRUE)))
    }
  }

  Xout
}

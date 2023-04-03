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

trapz <- function(x, y, dims = 1) {
  if ((dims - 1) > 0)
    perm <- c(dims:max(ndims(y), dims), 1:(dims - 1))
  else
    perm <- c(dims:max(ndims(y), dims))

  if (ndims(y) == 0)
    m <- 1
  else {
    if (length(x) != dim(y)[dims])
      cli::cli_abort("Dimensions mismatch in trapezoidal integration (function {.fn trapz}).")
    y = aperm(y, perm)
    m <- nrow(y)
  }

  if (m == 1) {
    M <- length(y)
    out <- sum(diff(x) * (y[-M] + y[-1]) / 2)
  } else {
    slice1 <- y[as.vector(outer(1:(m - 1), dim(y)[1] * (1:prod(dim(y)[-1]) - 1), '+'))]
    dim(slice1) <- c(m - 1, length(slice1) / (m - 1))
    slice2 <- y[as.vector(outer(2:m, dim(y)[1] * (1:prod(dim(y)[-1]) - 1), '+'))]
    dim(slice2) <- c(m - 1, length(slice2) / (m - 1))
    out <- t(diff(x)) %*% (slice1 + slice2) / 2.
    siz <- dim(y)
    siz[1] <- 1
    out <- array(out, siz)
    perm2 <- rep(0, length(perm))
    perm2[perm] <- 1:length(perm)
    out <- aperm(out, perm2)
    ind <- which(dim(out) != 1)
    out <- array(out, dim(out)[ind])
  }

  out
}

ndims <- function(x){
  length(dim(x))
}

check_centroid_type <- function(x) {
  if (x == "mean" || x == "medoid") return(list(name = x, extra = 0))
  if (grepl("lowess", x, fixed = TRUE)) {
    span_value <- as.numeric(gsub("\\D", "", x)) / 100
    if (is.na(span_value)) span_value <- 0.1
    return(list(name = "lowess", extra = span_value))
  }
  if (grepl("poly", x, fixed = TRUE)) {
    degree_value <- as.numeric(gsub("\\D", "", x))
    if (is.na(degree_value)) degree_value <- 4
    return(list(name = "lowess", extra = degree_value))
  }
  cli::cli_abort("The input argument {.arg centroid_type} should be one of {.code mean}, {.code medoid}, {.code lowessXX} or {.code polyXX}.")
}

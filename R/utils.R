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
  if (x == "mean" || x == "median" || x == "medoid") return(list(name = x, extra = 0))
  if (grepl("lowess", x, fixed = TRUE)) {
    span_value <- as.numeric(gsub("\\D", "", x)) / 100
    if (is.na(span_value)) span_value <- 0.1
    return(list(name = "lowess", extra = span_value))
  }
  if (grepl("poly", x, fixed = TRUE)) {
    degree_value <- as.numeric(gsub("\\D", "", x))
    if (is.na(degree_value)) degree_value <- 4
    return(list(name = "poly", extra = degree_value))
  }
  cli::cli_abort("The input argument {.arg centroid_type} should be one of {.code mean}, {.code medoid}, {.code lowessXX} or {.code polyXX}.")
}

format_inputs <- function(x, y = NULL) {
  # Here x is N x M and y is N x L x M when provided
  if (is.null(y) && rlang::is_installed("funData")) {
    if (inherits(x, "funData")) {
      if (length(x@argvals) > 1)
        cli::cli_abort(c(
          "The {.pkg fdacluster} package does not support functional data defined ",
          "on multivariate domains."
        ))
      L <- 1
      y <- x@X
      dims <- dim(y)
      N <- dims[1]
      M <- dims[2]
      y <- array(y, dim = c(N, L, M))
      x <- x@argvals[[1]]
    } else if (inherits(x, "irregFunData")) {
      L <- 1
      N <- length(x@argvals)
      M <- x@argvals |>
        purrr::map_int(length) |>
        mean() |>
        round()
      y <- array(dim = c(N, L, M))
      y[, 1, ] <- x@X |>
        purrr::imap(\(values, id) stats::approx(x@argvals[[id]], values, n = M)$y) |>
        do.call(rbind, args = _)
      x <- x@argvals |>
        purrr::map(\(grid) seq(min(grid), max(grid), length.out = M)) |>
        do.call(rbind, args = _)
    } else if (inherits(x, "multiFunData")) {
      L <- length(x)
      dims <- dim(x[[1]]@X)
      grid <- x[[1]]@argvals[[1]]
      purrr::walk(x, \(fData) {
        if (length(fData@argvals) != 1)
          cli::cli_abort(c(
            "The {.pkg fdacluster} package does not support functional data ",
            "defined on multivariate domains."
          ))
        if (any(fData@argvals[[1]] != grid))
          cli::cli_abort(c(
            "All components of the {.cls multiFunData} object must share the same ",
            "evaluation grids."
          ))
        if (any(dim(fData@X) != dims))
          cli::cli_abort(c(
            "All components of the {.cls multiFunData} object must have values ",
            "stored in matrices with the same dimensions."
          ))
      })
      N <- dims[1]
      M <- dims[2]
      y <- array(dim = c(N, L, M))
      for (l in 1:L) y[, l, ] <- x[[l]]@X
      x <- grid
    } else
      cli::cli_abort(c(
        "Functional data provided in a single argument {.arg x} must be either of ",
        "class {.cls funData} or of class {.cls irregFunData} or of class ",
        "{.cls multiFunData}."
      ))
  } else if (rlang::is_installed("fda") && inherits(y, "fd")) {
    dims <- purrr::map_int(y$fdnames, length)
    M <- dims[1]
    N <- dims[2]
    L <- dims[3]
    if (is.vector(x)) {
      if (length(x) != M)
        cli::cli_abort(c(
          "The number of function evaluations ({M}) does not match the grid ",
          "size ({length(x)})."
        ))
      y <- fda::eval.fd(x, y)
    } else {
      if (nrow(x) != N)
        cli::cli_abort(c(
          "When provided multiple evaluation grids as a matrix, the number of ",
          "rows should match the number of curves."
        ))
      if (ncol(x) != M)
        cli::cli_abort(c(
          "When provided multiple evaluation grids as a matrix, the number of ",
          "columns should match the common grid size."
        ))
      y <- fda::eval.fd(t(x), y)
    }
    if (is.null(dim(y))) {
      # y is a single 1-dimensional curve
      y <- array(y, dim = c(length(y), 1, 1))
    } else if (length(dim(y)) == 2) {
      # y is N 1-dimensional curves
      y <- array(y, dim = c(dim(y), 1))
    }
    y <- aperm(y, c(2, 3, 1))
  } else {
    if (length(dim(y)) == 2) {
      y <- array(y, c(dim(y)[1], 1, dim(y)[2]))
    }
    dims <- dim(y)
    N <- dims[1]
    L <- dims[2]
    M <- dims[3]
  }

  # Handle vector grid
  if (is.vector(x)) {
    x <- matrix(x, N, M, byrow = TRUE)
  }

  if (anyNA(x))
    cli::cli_abort("The input argument {.arg x} should not contain non-finite values.")

  if (anyNA(y))
    cli::cli_abort("The input argument {.arg y} should not contain non-finite values.")

  # output x matrix NxM and y array NxLxM
  list(x = x, y = y)
}

remove_missing_points <- function(grids, curves) {
  dims <- dim(curves)
  N <- dims[1]
  L <- dims[2]
  M <- dims[3]
  out_grids <- grids
  out_curves <- curves
  for (n in 1:N) {
    non_na_indices <- !is.na(curves[n, 1, ])
    tmin <- min(grids[n, non_na_indices])
    tmax <- max(grids[n, non_na_indices])
    tout <- seq(tmin, tmax, length.out = M)
    out_grids[n, ] <- tout
    for (l in 1:L)
      out_curves[n, l, ] <- stats::approx(grids[n, ], curves[n, l, ], xout = tout)$y
  }
  list(grids = out_grids, curves = out_curves)
}

prototype_curve <- function(s, a = 1, b = 1) {
  matrix(a, nrow = nrow(s), ncol = ncol(s)) * sin(s) +
    matrix(b, nrow = nrow(s), ncol = ncol(s)) * sin(s^2 / (2 * pi))
}

n_total <- 90
p <- 100
set.seed(1234)
error_terms <- matrix(rnorm(4 * n_total, sd = 0.05), nrow = n_total)
x_base <- seq(0, 2 * pi, length.out = p)
x_perturbed <- error_terms[, 3] + (1 + error_terms[, 4]) %*% matrix(x_base, nrow = 1)
x_deformed <- -1/3 + error_terms[, 3] + (3/4 + error_terms[, 4]) %*% matrix(x_base, nrow = 1)

png(filename = "data-raw/simulated90.png")
par(mfrow = c(2, 2))

# Case A ------------------------------------------------------------------

y_caseA <- prototype_curve(x_perturbed, 1 + error_terms[, 1], 1 + error_terms[, 1])
matplot(x_base, t(y_caseA), type = "l", main = "Case A", xlab = "s", ylab = "")

# Case B ------------------------------------------------------------------

y_caseB <- y_caseA
idx <- 46:90
y_caseB[idx, ] <- prototype_curve(x_perturbed[idx, ], 2 + error_terms[idx, 1], -1 + error_terms[idx, 1])
matplot(x_base, t(y_caseB), type = "l", main = "Case B", xlab = "s", ylab = "")

# Case C ------------------------------------------------------------------

y_caseC <- y_caseA
idx <- 46:90
y_caseC[idx, ] <- prototype_curve(x_deformed[idx, ], 1 + error_terms[idx, 1], 1 + error_terms[idx, 1])
matplot(x_base, t(y_caseC), type = "l", main = "Case C", xlab = "s", ylab = "")

# Case D ------------------------------------------------------------------

y_caseD <- y_caseA
idx <- 31:60
y_caseD[idx, ] <- y_caseC[idx, ]
idx <- 61:90
y_caseD[idx, ] <- y_caseB[idx, ]
matplot(x_base, t(y_caseD), type = "l", main = "Case D", xlab = "s", ylab = "")

dev.off()

# Data to be saved for use in kma
simulated90 <- list()
simulated90$x <- x_base
simulated90$y <- y_caseD

usethis::use_data(simulated90, overwrite = TRUE, compress = "xz")

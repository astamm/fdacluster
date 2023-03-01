library(magrittr)
library(tibble)

# Download data from webpage
tmp_dir <- tempdir()
zipfile <- paste0(tmp_dir, "/aneurisk65.zip")
zipdir <- paste0(tmp_dir, "/aneurisk")
download.file(
  url = "https://statistics.mox.polimi.it/wp-content/uploads/2015/07/AneuRisk65_ReadMe.zip",
  destfile = zipfile
)
unzip(
  zipfile = zipfile,
  exdir = zipdir
)

# Import all data in R
all_data <- tibble()
for (i in 1:65) {
  df <- readr::read_delim(
    file = paste0(zipdir, "/Rawdata_FKS_", i, ".txt"),
    delim = " "
  )
  all_data <- dplyr::bind_rows(
    all_data,
    df |> dplyr::mutate(Patient = i)
  )
}

unlink(tmp_dir, recursive = TRUE)

# Select only grids and first derivatives
selected_data <- all_data |>
  dplyr::select(
    Patient,
    s = Curv_Abscissa,
    x = X1_FKS_ref,
    y = Y1_FKS,
    z = Z1_FKS
  ) |>
  tidyr::nest(
    grids = s,
    values = c(x, y, z)
  )

# Put grids and values in array format
n <- nrow(selected_data)
d <- unique(purrr::map_int(selected_data$values, ncol))
p <- max(purrr::map_int(selected_data$grids, nrow))
aneurisk_grids <- matrix(nrow = n, ncol = p)
aneurisk_values <- array(dim = c(n, d, p))

for (i in 1:n) {
  s <- selected_data$grids[[i]]$s
  sout <- seq(
    from = min(s, na.rm = TRUE),
    to = max(s, na.rm = TRUE),
    length.out = p
  )
  x <- selected_data$values[[i]]$x
  y <- selected_data$values[[i]]$y
  z <- selected_data$values[[i]]$z
  x <- approx(s, x, sout)$y
  y <- approx(s, y, sout)$y
  z <- approx(s, z, sout)$y
  aneu_grids[i, ] <- sout
  aneu_values[i, , ] <- rbind(x, y, z)
}

aneurisk65 <- list()
aneurisk65$x <- aneurisk_grids
aneurisk65$y <- aneurisk_values

usethis::use_data(aneurisk65, overwrite = TRUE, compress = "xz")

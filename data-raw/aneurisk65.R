library(magrittr)
library(tibble)

# Download data from webpage
tmp_dir <- tempdir()
zipfile <- paste0(tmp_dir, "/aneurisk65.zip")
zipdir <- paste0(tmp_dir, "/aneurisk")
download.file(
  url = "http://ecm2.mathcs.emory.edu/aneuriskdata/files/Carotid-data_MBI_workshop.zip",
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

xValues <- funData::irregFunData(
  argvals = purrr::map(selected_data$grids, "s"),
  X = purrr::map(selected_data$values, "x")
)

out <- fdakmeans(xValues, n_clusters = 2)

funData::plot(yValues)

M <- xValues@argvals |> map_int(length) |> mean() |> round()
xValues@argvals |>
  map(\(grid) seq(min(grid), max(grid), length.out = M)) |>
  do.call(rbind, args = _) |> dim()
xValues@X |>
  imap(\(values, id) approx(xValues@argvals[[id]], values, n = M)$y) |>
  do.call(rbind, args = _) |> dim()

yValues <- funData::irregFunData(
  argvals = purrr::map(selected_data$grids, "s"),
  X = purrr::map(selected_data$values, "y")
)

zValues <- funData::irregFunData(
  argvals = purrr::map(selected_data$grids, "s"),
  X = purrr::map(selected_data$values, "z")
)

allValues <- funData::multiFunData(xValues, yValues, zValues)

growthData <- funData::funData(argvals = fda::growth$age, X = t(fda::growth$hgtm))
dim(growthData@X)

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
  x <- stats::approx(s, x, sout)$y
  y <- stats::approx(s, y, sout)$y
  z <- stats::approx(s, z, sout)$y
  aneu_grids[i, ] <- sout
  aneu_values[i, , ] <- rbind(x, y, z)
}

aneurisk65 <- list()
aneurisk65$x <- aneurisk_grids
aneurisk65$y <- aneurisk_values

# usethis::use_data(aneurisk65, overwrite = TRUE, compress = "xz")

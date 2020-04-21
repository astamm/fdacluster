
<!-- README.md is generated from README.Rmd. Please edit that file -->

# fdakmapp

<!-- badges: start -->

[![R build
status](https://github.com/astamm/fdakmapp/workflows/R-CMD-check/badge.svg)](https://github.com/astamm/fdakmapp/actions)
[![Coverage](https://github.com/astamm/fdakmapp/workflows/test-coverage/badge.svg)](https://github.com/astamm/fdakmapp/actions)
[![Website](https://github.com/astamm/fdakmapp/workflows/pkgdown/badge.svg)](https://github.com/astamm/fdakmapp/actions)
<!-- badges: end -->

The **fdakmapp** package provides the `kma()` function that jointly
performs clustering and alignment of a functional data set. Mean and
medoid center methods are available. The algorithm is fully
parallelized.

## Getting Started

The source code can be cloned or downloaded directly from github. An R
studio project file is provided to open the project in RStudio.

### Prerequisites

The package is linked against OpenMP, the BLAS and LAPACK libraries.
Therefore, you should make sure that your R / RStudio environment is
ready

  - to build packages: to check this, type in the console

<!-- end list -->

    # install.packages("pkgbuild")
    pkgbuild::has_devel()

  - to build packages with compiled code:
    
      - on Linux: it should be ready;
      - on Windows: you should install
        [Rtools](https://cloud.r-project.org/bin/windows/Rtools/);
      - on macOS: you should install a [set of
        tools](https://cloud.r-project.org/bin/macosx/tools/).

  - to build packages with C++ code through Rcpp: on macOS, you can for
    instance follow this
    [tutorial](https://thecoatlessprofessor.com/programming/cpp/r-compiler-tools-for-rcpp-on-macos/).

### Installation

The **fdakmapp** package is currently only available on GitHub. You can
install it as follows:

    # install.packages("remotes")
    remotes::install_github("astamm/fdakmapp")

### Example

``` r
library(fdakmapp)
res <- kma(
  simulated30$x,
  simulated30$y,
  seeds = c(21, 13),
  n_clust = 2,
  center_method = "medoid",
  warping_method = "affine",
  dissimilarity_method = "pearson"
)
#> Information about the data set:
#>  - Number of observations: 30
#>  - Number of dimensions: 1
#>  - Number of points: 200
#> Information about cluster initialization:
#>  - Number of clusters: 2
#>  - Initial seeds for cluster centers:         21        13
#> 
#> Information about the methods used within the algorithm:
#>  - Warping method: affine
#>  - Center method: medoid
#>  - Dissimilarity method: pearson
#>  - Optimization method: bobyqa
#> Information about warping parameter bounds:
#>  - Warping options:    0.1500   0.1500
#> 
#> Information about convergence criteria:
#>  - Maximum number of iterations: 100
#>  - Distance relative tolerance: 0.001
#> Information about parallelization setup:
#>  - Number of threads: 1
#>  - Parallel method: 0
#> Other information:
#>  - Use fence to robustify: 0
#>  - Check total dissimilarity: 1
#>  - Compute overall center: 0
#> Start execution.
#> Compute initial templates: Done.
#> Start while iteration
#> Iteration #1
#> Set bound: Done
#> 1. Compute best warping: 
#> Done
#> Cluster indices: 
#>         0        1
#>  - Size of cluster #0: 10
#>  - Size of cluster #1: 20
#> Parameter normalization: Done
#> Update individual warped grids: Done
#> Compute new templates: 
#> Check total similarity: Done
#> Iteration #2
#> Set bound: Done
#> 2. Compute best warping: 
#> Done
#> Cluster indices: 
#>         0        1
#>  - Size of cluster #0: 10
#>  - Size of cluster #1: 20
#> Parameter normalization: Done
#> Update individual warped grids: Done
#> Compute new templates: 
#> Check total similarity: Done
#> End while iterations
#> Final warping: Done

tmp1 <- res$y
tmp2 <- res$y_centers_final
res$y <- array(dim = c(dim(tmp1)[1], dim(tmp1)[3], dim(tmp1)[2]))
res$y_centers_final <- array(dim = c(dim(tmp2)[1], dim(tmp2)[3], dim(tmp2)[2]))
for (i in 1) {
  res$y[, , i] <- tmp1[, i, ]
  res$y_centers_final[, , i] <- tmp2[, i, ]
}

plot(res)
```

<img src="man/figures/README-example-simulated30-1.png" width="100%" />

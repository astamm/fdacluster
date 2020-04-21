<!-- badges: start -->
[![R build status](https://github.com/astamm/fdakmapp/workflows/R-CMD-check/badge.svg)](https://github.com/astamm/fdakmapp/actions)
[![Coverage](https://github.com/astamm/fdakmapp/workflows/test-coverage/badge.svg)](https://github.com/astamm/fdakmapp/actions)
[![Website](https://github.com/astamm/fdakmapp/workflows/pkgdown/badge.svg)](https://github.com/astamm/fdakmapp/actions)
<!-- badges: end -->

# fdakmapp
The fdakmapp package provides the kmap function that jointly performs clustering and alignment of a functional 
dataset (multidimensional or unidimensional). The centers can be computed by mean or medoid center methods.
The parallel version is available.

## Getting Started
The source code can be cloned or downloaded directly from github.
An R studio project file is provided to open the project in RStudio.

### Prerequisites

The package is linked against OpenMP, the BLAS and LAPACK libraries in Makevars.

To install the package locally the packages Rcpp and RcppArmadillo needs to be installed.
Using install_github() will install the dependencies (Rcpp and RcppArmadillo) automatically.

### Installing

The package can be installed directly from github but devtools is required.

If devtools is not installed use the following comand to install it
```
install.packages('devtools') 
```
and then install the package
```
library(devtools)
install_github('zitale/fdakmapp')
```

Usually the clag compiler doesn't support OpenMP. For this reason is available a release without -fopenmp.
```
library(devtools)
install_github('zitale/fdakmapp@v2.0.2.noomp')
```

Othewise the problem can be solved following the tutorial at: https://clang-omp.github.io/.

### Automatic tests

```
devtools::test()
```
### Example

```
res<-kmap(x=aneurisk65$x, y=aneurisk65$y, n_clust=2)
kmap_show_results(res,FALSE,FALSE)
```

## Documentation

The R documentation can be found in the main directory in fdakmapp.pdf.
The C++ documentation can be found at https://zitale.github.io/fdakmapp/.

# fdakmapp
The fdakmapp package provides the kmap function that jointly performs clustering and alignment of a functional 
dataset (multidimensional or unidimensional). The centers can be computed by mean and medoid center methods.
Many options are available and also the parallal version.

## Getting Started
The package can be cloned or downloaded directly from github.
An R studio project file is provided to open the project in RStudio.

### Prerequisites

The package is linked against OpenMP, the BLAS and LAPACK libraries, to use RcppArmadillo, in Makevars.

To install the package locally the packages Rcpp and RcppArmadillo needs to be installed.
Using install_github() will install the dependencies (Rcpp and RcppArmadillo) automatically.

### Installing

The package can be installed directly from github but devtools is required.

```
library(devtools)
install_github('zitale/fdakmapp')
```
Usally the Mac's compiler doesn't support OpenMP. For this reason is available a release the compile the package without -fopenmp.

```
library(devtools)
install_github('zitale/fdakmapp@v2.0.2.noomp')
```


Othewise another compiler can be installed following the tutorial at: https://clang-omp.github.io/.

### Automatic tests

```
devtools::test()
```
### Example

```
res<-kmap(x=aneurisk65$x, y=aneurisk65$y, n_clust=2)
kmap_show_results(res,FALSE)
```

## Documentation

The R documentation can be found in the main directory in fdakmapp.pdf.
The C++ documentation can be found at https://zitale.github.io/fdakmapp/.

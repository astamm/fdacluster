# fdakmapp
The fdakmapp package provides the kmap function that jointly performs clustering and alignment of a functional 
dataset (multidimensional or unidimensional). The centers can be computed by mean and medoid center methods.
Many options are available and also the parallal version.

## Getting Started
The package can be cloned and downloaded directly from github.
An R studio project file is provided to open the project in RStudio.

### Prerequisites

To install the package locally the packages Rcpp and RcppArmadillo needs to be instaleld.
Using install_github() will install the dependencies automatically.

### Installing

The package can be installed directly from github but devtools is required.

```
install.packages(devtools)
library(devtools)
install_github('zitale/fdakmapp')
```

### Automatic tests

```
devtools::test()
```
### Example

```
res<-kmap(x=aneurisk65$x, y=aneurisk65$y, n_clust=2)
kmap_show_results(res,FALSE)
```

# An `mcaps` object from simulated data for examples

An object of class `mcaps` storing the result of the
[`compare_caps()`](https://astamm.github.io/fdacluster/reference/compare_caps.md)
function applied on the data set
[`simulated30_sub`](https://astamm.github.io/fdacluster/reference/simulated30_sub.md)
for comparing the clustering structures found by the
[`fdakmeans()`](https://astamm.github.io/fdacluster/reference/fdakmeans.md)
function with `mean` centroid type used with various classes of warping
functions and varying number of clusters.

## Usage

``` r
sim30_mcaps
```

## Format

An object of class `mcaps` which is effectively a
[tibble::tibble](https://tibble.tidyverse.org/reference/tibble.html)
with 5 columns and as many rows as there are clustering strategies to
compare. The 5 column-variables are:

- `n_clusters`: The number of clusters;

- `clustering_method`: The clustering method;

- `warping_class`: The class of warping functions used for curve
  alignment;

- `centroid_type`: The type of centroid used to compute a cluster
  representative;

- `caps_obj`: The result of the corresponding clustering strategy as
  objects of class
  [`caps`](https://astamm.github.io/fdacluster/reference/caps.md).

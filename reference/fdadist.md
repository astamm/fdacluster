# Computes the distance matrix for functional data with amplitude and phase separation

This function computes the matrix of pairwise distances between curves a
functional data sample. This can be achieved with or without phase and
amplitude separation, which can be done using a variety of warping
classes.

## Usage

``` r
fdadist(
  x,
  y = NULL,
  is_domain_interval = FALSE,
  transformation = c("identity", "srvf"),
  warping_class = c("none", "shift", "dilation", "affine", "bpd"),
  metric = c("l2", "normalized_l2", "pearson"),
  cluster_on_phase = FALSE,
  labels = NULL
)
```

## Arguments

- x:

  A numeric vector of length \\M\\ or a numeric matrix of shape \\N
  \times M\\ or an object of class
  [`funData::funData`](https://rdrr.io/pkg/funData/man/funData-class.html).
  If a numeric vector or matrix, it specifies the grid(s) of size \\M\\
  on which each of the \\N\\ curves have been observed. If an object of
  class
  [`funData::funData`](https://rdrr.io/pkg/funData/man/funData-class.html),
  it contains the whole functional data set and the `y` argument is not
  used.

- y:

  Either a numeric matrix of shape \\N \times M\\ or a numeric array of
  shape \\N \times L \times M\\ or an object of class
  [`fda::fd`](https://rdrr.io/pkg/fda/man/fd.html). If a numeric matrix
  or array, it specifies the \\N\\-sample of \\L\\-dimensional curves
  observed on grids of size \\M\\. If an object of class
  [`fda::fd`](https://rdrr.io/pkg/fda/man/fd.html), it contains all the
  necessary information about the functional data set to be able to
  evaluate it on user-defined grids.

- is_domain_interval:

  A boolean specifying whether the sample of curves is defined on a
  fixed interval. Defaults to `FALSE`.

- transformation:

  A string specifying the transformation to apply to the original sample
  of curves. Choices are no transformation
  (`transformation = "identity"`) or square-root velocity function
  `transformation = "srvf"`. Defaults to `"identity"`.

- warping_class:

  A string specifying the class of warping functions. Choices are no
  warping (`warping_class = "none"`), shift `y = x + b`
  (`warping_class = "shift"`), dilation `y = ax`
  (`warping_class = "dilation"`), affine `y = ax + b`
  (`warping_class = "affine"`) or boundary-preserving diffeomorphism
  (`warping_class = "bpd"`). Defaults to `"none"`.

- metric:

  A string specifying the metric used to compare curves. Choices are
  `"l2"`, `"normalized_l2"` or `"pearson"`. If
  `transformation == "srvf"`, the metric **must be** `"l2"` because the
  SRVF transform maps absolutely continuous functions to
  square-integrable functions. If `transformation == "identity"` and
  `warping_class` is either `dilation` or `affine`, the metric cab be
  either `"normalized_l2"` or `"pearson"`. The L2 distance is indeed
  **not** dilation-invariant or affine-invariant. The metric can also be
  `"l2"` if `warping_class == "shift"`. Defaults to `"l2"`.

- cluster_on_phase:

  A boolean specifying whether clustering should be based on phase
  variation or amplitude variation. Defaults to `FALSE` which implies
  amplitude variation.

- labels:

  A character vector specifying curve labels. Defaults to `NULL` which
  uses sequential numbers as labels.

## Value

A [stats::dist](https://rdrr.io/r/stats/dist.html) object storing the
distance matrix between the input curves using the metric specified
through the argument `metric` and the warping class specified by the
argument `warping_class`.

## Examples

``` r
idx <- c(1:5, 11:15, 21:25)
D <- fdadist(simulated30_sub$x[idx, ], simulated30_sub$y[idx, , ])
```

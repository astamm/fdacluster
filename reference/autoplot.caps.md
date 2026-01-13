# Visualizes the result of a clustering strategy stored in a `caps` object with ggplot2

This function creates a visualization of the result of the k-mean
alignment algorithm and invisibly returns the corresponding
[ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object which enable further customization of the plot. The user can
choose to visualize either the amplitude information data in which case
original and aligned curves are shown or the phase information data in
which case the estimated warping functions are shown.

## Usage

``` r
# S3 method for class 'caps'
autoplot(object, type = c("amplitude", "phase"), ...)
```

## Arguments

- object:

  An object of class
  [`caps`](https://astamm.github.io/fdacluster/reference/caps.md).

- type:

  A string specifying the type of information to display. Choices are
  `"amplitude"` for plotting the original and aligned curves which
  represent amplitude information data or `"phase"` for plotting the
  corresponding warping functions which represent phase information
  data. Defaults to `"amplitude"`.

- ...:

  Not used.

## Value

A [ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html)
object invisibly.

## Examples

``` r
ggplot2::autoplot(sim30_caps, type = "amplitude")

ggplot2::autoplot(sim30_caps, type = "phase")
```

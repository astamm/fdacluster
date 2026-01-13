# Diagnostic plot for the result of a clustering strategy stored in a `caps` object

This function plots the values of the distance to center and silhouette
for each observation. Observations are ordered within cluster by
decreasing value of silhouette.

## Usage

``` r
diagnostic_plot(x)
```

## Arguments

- x:

  An object of class
  [`caps`](https://astamm.github.io/fdacluster/reference/caps.md).

## Value

An object of class
[ggplot2::ggplot](https://ggplot2.tidyverse.org/reference/ggplot.html).

## Examples

``` r
diagnostic_plot(sim30_caps)
```

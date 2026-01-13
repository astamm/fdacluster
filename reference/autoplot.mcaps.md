# Visualizes results of multiple clustering strategies using ggplot2

This is an S3 method implementation of the
[`ggplot2::autoplot()`](https://ggplot2.tidyverse.org/reference/autoplot.html)
generic for objects of class `mcaps` to visualize the performances of
multiple [`caps`](https://astamm.github.io/fdacluster/reference/caps.md)
objects applied on the same data sets either in terms of WSS or in terms
of silhouette values.

## Usage

``` r
# S3 method for class 'mcaps'
autoplot(
  object,
  validation_criterion = c("wss", "silhouette"),
  what = c("mean", "distribution"),
  ...
)
```

## Arguments

- object:

  An object of class `mcaps`.

- validation_criterion:

  A string specifying the validation criterion to be used for the
  comparison. Choices are `"wss"` or `"silhouette"`. Defaults to
  `"wss"`.

- what:

  A string specifying the kind of information to display about the
  validation criterion. Choices are `"mean"` (which plots the mean
  values) or `"distribution"` (which plots the boxplots). Defaults to
  `"mean"`.

- ...:

  Other arguments passed to specific methods.

## Value

An object of class
[`ggplot2::ggplot`](https://ggplot2.tidyverse.org/reference/ggplot.html).

## Examples

``` r
p <- ggplot2::autoplot(sim30_mcaps)
```

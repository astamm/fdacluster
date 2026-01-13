# Plots the result of a clustering strategy stored in a `caps` object

This function creates a visualization of the result of the k-mean
alignment algorithm **without** returning the plot data as an object.
The user can choose to visualize either the amplitude information data
in which case original and aligned curves are shown or the phase
information data in which case the estimated warping functions are
shown.

## Usage

``` r
# S3 method for class 'caps'
plot(x, type = c("amplitude", "phase"), ...)
```

## Arguments

- x:

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

## Examples

``` r
plot(sim30_caps, type = "amplitude")

plot(sim30_caps, type = "phase")
```

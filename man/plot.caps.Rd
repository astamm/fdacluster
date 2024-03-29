% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/caps-plot.R
\name{plot.caps}
\alias{plot.caps}
\title{Plots the result of a clustering strategy stored in a \code{caps} object}
\usage{
\method{plot}{caps}(x, type = c("amplitude", "phase"), ...)
}
\arguments{
\item{x}{An object of class \code{\link{caps}}.}

\item{type}{A string specifying the type of information to display. Choices
are \code{"amplitude"} for plotting the original and aligned curves which
represent amplitude information data or \code{"phase"} for plotting the
corresponding warping functions which represent phase information data.
Defaults to \code{"amplitude"}.}

\item{...}{Not used.}
}
\description{
This function creates a visualization of the result of the k-mean alignment
algorithm \strong{without} returning the plot data as an object. The user can
choose to visualize either the amplitude information data in which case
original and aligned curves are shown or the phase information data in which
case the estimated warping functions are shown.
}
\examples{
plot(sim30_caps, type = "amplitude")
plot(sim30_caps, type = "phase")
}

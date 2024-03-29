% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/caps-plot.R
\name{autoplot.caps}
\alias{autoplot.caps}
\title{Visualizes the result of a clustering strategy stored in a \code{caps} object with
ggplot2}
\usage{
\method{autoplot}{caps}(object, type = c("amplitude", "phase"), ...)
}
\arguments{
\item{object}{An object of class \code{\link{caps}}.}

\item{type}{A string specifying the type of information to display. Choices
are \code{"amplitude"} for plotting the original and aligned curves which
represent amplitude information data or \code{"phase"} for plotting the
corresponding warping functions which represent phase information data.
Defaults to \code{"amplitude"}.}

\item{...}{Not used.}
}
\value{
A \link[ggplot2:ggplot]{ggplot2::ggplot} object invisibly.
}
\description{
This function creates a visualization of the result of the k-mean alignment
algorithm and invisibly returns the corresponding \link[ggplot2:ggplot]{ggplot2::ggplot} object
which enable further customization of the plot. The user can choose to
visualize either the amplitude information data in which case original and
aligned curves are shown or the phase information data in which case the
estimated warping functions are shown.
}
\examples{
\dontshow{if (requireNamespace("ggplot2", quietly = TRUE)) (if (getRversion() >= "3.4") withAutoprint else force)(\{ # examplesIf}
ggplot2::autoplot(sim30_caps, type = "amplitude")
ggplot2::autoplot(sim30_caps, type = "phase")
\dontshow{\}) # examplesIf}
}

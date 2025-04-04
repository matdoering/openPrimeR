% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots_constraints.R
\name{plot_constraint.histogram}
\alias{plot_constraint.histogram}
\title{Histogram of Constraints.}
\usage{
plot_constraint.histogram(
  primer.df,
  con.cols,
  con.identifier,
  boundaries = NULL,
  x.limits = NULL
)
}
\arguments{
\item{primer.df}{Primer data frame, not necessarily a \code{Primers} object.}

\item{con.cols}{Constraint identifiers in \code{primer.df} to plot.}

\item{con.identifier}{Name of the constraint to plot.}

\item{boundaries}{List with constraint settings.}

\item{x.limits}{Interval limiting the extent of the x-axis.}
}
\value{
A constraint histogram plot.
}
\description{
Plots a histogram of constraint values.
}
\keyword{internal}

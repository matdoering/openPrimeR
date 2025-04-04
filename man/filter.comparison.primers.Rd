% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.R
\name{filter.comparison.primers}
\alias{filter.comparison.primers}
\title{Filter Multiple Primer Sets.}
\usage{
filter.comparison.primers(
  primers,
  templates,
  active.constraints,
  settings,
  updateProgress = NULL
)
}
\arguments{
\item{primers}{List with primer data frames.}

\item{templates}{List with template data frames.}

\item{active.constraints}{Strings giving the constraints that are to be checked.}

\item{settings}{List with settings.}

\item{updateProgress}{Progress callback function for shiny.}
}
\value{
A list with filtered primer data frames.
}
\description{
Filters multiple primer sets at once.
}
\keyword{internal}

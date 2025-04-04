% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization_ILP.R
\name{get.sets.from.decisions}
\alias{get.sets.from.decisions}
\title{Optimal Sets from Decision Variables}
\usage{
get.sets.from.decisions(ILP.df, Tm.sets)
}
\arguments{
\item{ILP.df}{Data frame with ILP optimization results.}

\item{Tm.sets}{List with primer data frames for every target melting temperature.}
}
\value{
A list with optimal primer data sets for every target temperature.
}
\description{
Determines primer sets from decision variables from ILP.
}
\keyword{internal}

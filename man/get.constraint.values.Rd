% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/filters.R
\name{get.constraint.values}
\alias{get.constraint.values}
\title{Get the Values of a Constraint.}
\usage{
get.constraint.values(con.name, cur.candidates, mode.directionality)
}
\arguments{
\item{con.name}{The name of the constraint.}

\item{cur.candidates}{The \code{Primers} data frame where the values should be retrieved.}

\item{mode.directionality}{The direction for which values should be retrieved.}
}
\value{
The constraint values corresponding to \code{con.name} for the primers \code{cur.candidates}.
}
\description{
Get the Values of a Constraint.
}
\keyword{internal}

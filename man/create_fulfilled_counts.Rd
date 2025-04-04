% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primer_significance.R
\name{create_fulfilled_counts}
\alias{create_fulfilled_counts}
\title{Creation of Fulfilled/Failed Constraint Counts.}
\usage{
create_fulfilled_counts(primer.df, eval.cols = NULL)
}
\arguments{
\item{primer.df}{An evaluated \code{Primers} object.}

\item{eval.cols}{Evaluation columns in \code{primer.df} to consider.
By default (\code{NULL}) all evaluation columns are considered.}
}
\value{
A data frame with the number of fulfilled/failed constraints
for \code{primer.df}.
}
\description{
Creates counts of fullfilled/failed constraints.
}
\keyword{internal}

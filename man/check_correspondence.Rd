% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization_greedy.R
\name{check_correspondence}
\alias{check_correspondence}
\title{Check of Primer and Template Correspondence.}
\usage{
check_correspondence(primer.df, template.df)
}
\arguments{
\item{primer.df}{An object of class \code{Primers}.}

\item{template.df}{An object of class \code{Templates}.}
}
\value{
\code{TRUE} if the primers and templates seem to correspond,
\code{FALSE} otherwise.
}
\description{
Checks whether the primers relate to the correct templates.
}
\keyword{internal}

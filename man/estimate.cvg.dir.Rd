% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primers.R
\name{estimate.cvg.dir}
\alias{estimate.cvg.dir}
\title{Estimation of Primer Coverage.}
\usage{
estimate.cvg.dir(seqs, k, id = "")
}
\arguments{
\item{seqs}{A character vector of sequences to evaluate coverage for.}

\item{k}{A numeric giving the size of the primers.}

\item{id}{An optional identifier for the primers.}
}
\value{
A data frame with binding information.
}
\description{
Estimates the possible coverage of primers using
probes of size \code{k} and only considering perfect
matches without consideration of ambiguities.
}
\keyword{internal}

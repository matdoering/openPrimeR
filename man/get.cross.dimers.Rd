% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/con_dimerization.R
\name{get.cross.dimers}
\alias{get.cross.dimers}
\title{Cross dimers}
\usage{
get.cross.dimers(
  primers.1,
  primers.2,
  ions,
  annealing.temp,
  check.idx = NULL,
  no.structures = FALSE,
  mode = c("symmetric", "asymmetric")
)
}
\arguments{
\item{primers.1}{Input primers.}

\item{primers.2}{Input primers.}

\item{ions}{Sodium-equivalent ionic concentration.}

\item{annealing.temp}{The PCR annealing temperature.}

\item{check.idx}{indices of primers for checking cross-dimerization}

\item{no.structures}{Whether to compute structures of dimers.}

\item{mode}{'symmetric', if \code{primers.1} and \code{primers.2}
carry the same information (i.e. fw-fw, rev-rev, fw-rev), 'asymetric' else.}
}
\value{
Data frame with potential cross dimers.
}
\description{
Computes all possible primer cross-dimers.
}
\keyword{internal}

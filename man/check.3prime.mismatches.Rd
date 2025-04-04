% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/con_primer_coverage.R
\name{check.3prime.mismatches}
\alias{check.3prime.mismatches}
\title{3' Mismatch Check.}
\usage{
check.3prime.mismatches(
  template.df,
  primer.df,
  mode.directionality = c("fw", "rev", "both")
)
}
\arguments{
\item{template.df}{Template data frame.}

\item{primer.df}{Primer data frame.}

\item{mode.directionality}{Primer directionality.}
}
\value{
Returns the distance of mismatches from the 3' terminal end of primers.
}
\description{
Check for mismatches at primer 3' ends.
}
\keyword{internal}

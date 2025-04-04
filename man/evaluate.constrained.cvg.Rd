% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/con_primer_coverage.R
\name{evaluate.constrained.cvg}
\alias{evaluate.constrained.cvg}
\title{Evaluation of Primer Coverage.}
\usage{
evaluate.constrained.cvg(
  template.df,
  primer.df,
  cvg.df,
  mode.directionality = c("fw", "rev", "both"),
  settings,
  updateProgress = NULL
)
}
\arguments{
\item{template.df}{Template data frame.}

\item{primer.df}{Primer data frame.}

\item{cvg.df}{Data frame with basic coverage entries.}

\item{mode.directionality}{Primer directionality.}

\item{settings}{A \code{DesignSettings} object.}

\item{updateProgress}{Progress callback function for shiny.}
}
\value{
Primer data frame with information on the covered template sequences.
}
\description{
Evaluates the coverage of a set of primers.
}
\keyword{internal}

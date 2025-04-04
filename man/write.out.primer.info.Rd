% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization_algo.R
\name{write.out.primer.info}
\alias{write.out.primer.info}
\title{Write Out Optimization Data}
\usage{
write.out.primer.info(
  opti.results.loc,
  optimal.primers.data,
  mode.directionality,
  settings,
  sample.name,
  template.df,
  max.degen
)
}
\arguments{
\item{opti.results.loc}{Folder where optimization data reside.}

\item{optimal.primers.data}{List with optimization results.}

\item{mode.directionality}{Direction of primers.}

\item{settings}{Settings used in the optimization procedure.
List containing fw, rev settings.}

\item{sample.name}{Name of template sample.}

\item{template.df}{Template data frame.}

\item{max.degen}{Maximal degeneracy of primers.}
}
\value{
Write-out of primer information to \code{opti.results.loc}.
}
\description{
Writes out all data relating to the optimization of primers.
}
\keyword{internal}

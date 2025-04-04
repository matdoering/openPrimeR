% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization_greedy.R
\name{select.best.primer.idx}
\alias{select.best.primer.idx}
\title{Greedy Choice}
\usage{
select.best.primer.idx(
  result,
  primers,
  deltaG.cutoff,
  target.temp,
  primer_conc,
  na_salt_conc,
  mg_salt_conc,
  k_salt_conc,
  tris_salt_conc
)
}
\arguments{
\item{result}{Data frame of current optimized primer data set that is to be augmented.}

\item{primers}{Data frame of candidate primers for addition to \code{result}.}

\item{deltaG.cutoff}{Free energy cutoff for cross-dimerization.}

\item{target.temp}{Target annealing temperature in Celsius.}

\item{primer_conc}{Primer concentration.}

\item{na_salt_conc}{Sodium ion concentration.}

\item{mg_salt_conc}{Magensium ion concentration.}

\item{k_salt_conc}{Potassium ion concentration.}

\item{tris_salt_conc}{Tris ion concentration.}
}
\value{
The index of a suitable primer according to Greedy selection.
}
\description{
Selects the currently best primer for Greedy primer selection.
}
\keyword{internal}

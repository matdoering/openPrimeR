% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/con_primer_efficiency.R
\name{compute.efficiency}
\alias{compute.efficiency}
\title{Primer Efficiency.}
\usage{
compute.efficiency(
  fw.primers,
  fw.start,
  fw.end,
  covered,
  taqEfficiency,
  annealing.temp,
  primer_conc,
  sodium.eq.concentration,
  mode.directionality = c("fw", "rev"),
  seqs
)
}
\arguments{
\item{fw.primers}{Primer sequence strings.}

\item{fw.start}{Binding position (start).}

\item{fw.end}{Binding position (end).}

\item{covered}{List of covered template indices per primer.}

\item{taqEfficiency}{Whether the efficiency shall be computed
using a mismatch-model developed for Taq polymerases. The default setting
is \code{TRUE}. Set \code{taqEfficiency} to \code{FALSE} if you are using
another polymerase than Taq.}

\item{annealing.temp}{Annealing temperature for which to evaluate efficiency.}

\item{primer_conc}{Primer concentration.}

\item{sodium.eq.concentration}{The sodium-equivalent concentration of ions.}

\item{mode.directionality}{Primer directionality.}

\item{seqs}{Template sequence strings.}
}
\value{
The efficiencies of primer binding events.
}
\description{
Computes the efficiency of primer binding events for Taq polymerase.
}
\details{
This function uses DECIPHER's \link[DECIPHER]{CalculateEfficiencyPCR}.
}
\references{
Wright, Erik S., et al. 
"Exploiting extension bias in polymerase chain reaction to improve 
primer specificity in ensembles of nearly identical DNA templates." 
Environmental microbiology 16.5 (2014): 1354-1365.
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/con_melting_temperature.R
\name{call.melt.single}
\alias{call.melt.single}
\title{Thermodynamic melting temperature computations.}
\usage{
call.melt.single(
  primers,
  complements,
  out.file,
  primer_conc,
  na_salt_conc,
  mg_salt_conc,
  k_salt_conc,
  tris_salt_conc,
  ID
)
}
\arguments{
\item{primers}{List of primer strings.}

\item{complements}{List with corresponding complements.}

\item{out.file}{Path to the file where MELTING will write the results.}

\item{primer_conc}{Primer concentration.}

\item{na_salt_conc}{Sodium ion concentration.}

\item{mg_salt_conc}{Magensium ion concentration.}

\item{k_salt_conc}{Potassium ion concentration.}

\item{tris_salt_conc}{Tris buffer concentration.}

\item{ID}{identifiers of the input primers}
}
\value{
Melting temperature data frame.
}
\description{
Computes the melting temperature for the input primers.
}
\keyword{internal}

% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization_ILP.R
\name{build.ILP.df}
\alias{build.ILP.df}
\title{Construction of ILP Results.}
\usage{
build.ILP.df(
  ILP,
  vars,
  primer.df,
  template.df,
  i,
  target.temp,
  time = NA,
  deltaG_Cutoff = NA,
  deltaG_Limit = NA
)
}
\arguments{
\item{ILP}{A solved ILP instance.}

\item{vars}{The ILP decision variables.}

\item{primer.df}{The primer data frame correspdong to the \code{ILP}.}

\item{template.df}{The template data frame.}

\item{i}{Index for the ILP.}

\item{target.temp}{Target melting temperature in Celsius.}

\item{time}{Runtime of the ILP.}

\item{deltaG_Cutoff}{Free energy cutoff used for the dimerization constraint.}

\item{deltaG_Limit}{The free energy boundary for dimerization.}
}
\value{
Data frame summarizing the ILP solution.
}
\description{
Constructs a data frame summarizing the properties of an ILP solution.
}
\keyword{internal}

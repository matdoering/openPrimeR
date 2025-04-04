% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/templates.R
\name{get.leader.exon.regions.single}
\alias{get.leader.exon.regions.single}
\title{Individual Binding Annotation}
\usage{
get.leader.exon.regions.single(
  l.seq,
  lex.seq,
  direction = c("fw", "rev"),
  gap.char
)
}
\arguments{
\item{l.seq}{Data frame with individual binding regions.}

\item{lex.seq}{Template data frame.}

\item{direction}{The primer direction for which the binding info is valid.}

\item{gap.char}{The character for gaps in alignments.}
}
\value{
Template data  frame with annotated binding regions.
}
\description{
Annotate individual binding regions.
}
\keyword{internal}

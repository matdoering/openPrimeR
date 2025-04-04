% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primers.R
\name{restriction_ali}
\alias{restriction_ali}
\title{Identification of Badly Fitting Regions.}
\usage{
restriction_ali(primer.seqs, template.seqs, search.hits)
}
\arguments{
\item{primer.seqs}{Primer sequences.}

\item{template.seqs}{Template sequences.}

\item{search.hits}{Template substrings that agree well
with the input primers.}
}
\value{
A list with putative restriction sites for every primer.
}
\description{
Identify regions in the templates where the primers are
not very complementary. These regions indicate possible
restriction enzyme adapters.
}
\keyword{internal}

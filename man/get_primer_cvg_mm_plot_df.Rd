% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/plots_coverage.R
\name{get_primer_cvg_mm_plot_df}
\alias{get_primer_cvg_mm_plot_df}
\title{Data for Mismatch Primer Coverage Plot.}
\usage{
get_primer_cvg_mm_plot_df(primer.df, template.df)
}
\arguments{
\item{primer.df}{A \code{Primers} object.}

\item{template.df}{A \code{Templates} object.}
}
\value{
A data frame for plotting mismatch primer coverage.
}
\description{
Ensures that there's an entry for every possible mismatch setting.
}
\keyword{internal}

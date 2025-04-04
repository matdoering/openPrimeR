% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/primers.R
\name{update_primer_binding_regions}
\alias{update_primer_binding_regions}
\title{Update of Primer Binding Regions.}
\usage{
update_primer_binding_regions(primer.df, template.df, old.template.df)
}
\arguments{
\item{primer.df}{A \code{Primers} data frame.}

\item{template.df}{Templates with the new binding regions.}

\item{old.template.df}{Templates with the old binding regions.}
}
\value{
A \code{Primers} object with updated relative binding positions.
}
\description{
Updates the relative primer binding sites in the templates
when the template binding regions have changed since the last
coverage computation.
}
\keyword{internal}

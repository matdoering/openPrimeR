% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/templates.R
\name{rbind.Templates}
\alias{rbind.Templates}
\title{rbind for Template class.}
\usage{
\method{rbind}{Templates}(...)
}
\arguments{
\item{...}{Parameters for rbind function.}
}
\value{
Row-binded Templates data frame.
}
\description{
Ensures that the rbind result has the appropriate class.
}
\examples{
data(Ippolito)
template.df <- rbind(template.df, template.df)
}
\keyword{internal}

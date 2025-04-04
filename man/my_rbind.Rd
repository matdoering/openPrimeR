% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helper_functions.R
\name{my_rbind}
\alias{my_rbind}
\title{Smartbind preserving classes.}
\usage{
my_rbind(...)
}
\arguments{
\item{...}{Data frames.}
}
\value{
A data frame resulting from row binding of \code{...}.
}
\description{
Rbind allowing for column mismatch, retains the classes of the data frames.
Motivation: smartbind/rbind.fill only keep the data.frame class
but not additional classes.
}
\keyword{internal}

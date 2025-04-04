% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/errors.R
\name{my.error}
\alias{my.error}
\title{Custom Error}
\usage{
my.error(subclass, message, call = sys.call(-1), ...)
}
\arguments{
\item{subclass}{String giving the specific type of error.}

\item{message}{Message to be displayed to the user.}

\item{call}{Environment where the error ocurred.}

\item{...}{Other arguments to be passed to  the condition function.}
}
\value{
Generates a custom error.
}
\description{
Creates an error with a custom class.
}
\keyword{internal}

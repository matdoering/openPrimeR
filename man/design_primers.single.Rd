% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/optimization_algo.R
\name{design_primers.single}
\alias{design_primers.single}
\title{Design Primers for a Single Direction}
\usage{
design_primers.single(
  template.df,
  sample.name,
  mode.directionality = c("fw", "rev"),
  settings,
  timeout,
  opti.algo,
  allowed.region.definition,
  init.algo,
  max.degen,
  conservation,
  target.temps,
  required.cvg,
  fw.primers = NULL,
  cur.results.loc = NULL,
  primer.df = NULL,
  updateProgress = NULL
)
}
\arguments{
\item{template.df}{Template data frame with sequences for which primers shall be designed.}

\item{sample.name}{Identifier for the templates.}

\item{mode.directionality}{Template strands for which primers shall be designed.
Primers can be designed either only for forward strands, only for reverse strands, or
both strand directions.}

\item{settings}{A \code{DesignSettings} object specifying the 
criteria for designing primers.}

\item{timeout}{Timeout in seconds for the optimization with ILPs.}

\item{opti.algo}{The algorithm to be used for solving the primer set covering problem.}

\item{allowed.region.definition}{Definition of the target binding sites used for evaluating the coverage.
If \code{allowed.region.definition} is "within", primers have to lie within the allowed binding region.
If \code{allowed.region.definition} is "any", primers have to overlap with the allowed binding region.
The default is that primers have to bind within the target binding region.}

\item{init.algo}{The algorithm to be used for initializing primers.
If \code{init.algo} is \code{naive}, then primers are constructed from substrings of the input template sequences.
If \code{init.algo} is \code{tree}, phylogenetic trees are used to form degenerate primers whose degeneracy is
bounded by \code{max.degen}.}

\item{max.degen}{Maximal degeneracy of merged primers.}

\item{conservation}{When using the tree-based primer initialization, consider only the \code{conservation}
percentile of regions with the highest conservation.}

\item{target.temps}{Target melting temperatures for optimized primer sets in Celsius.
Only required when optimizing primers for both strand directions and one optimization was already performed.}

\item{required.cvg}{The target ratio of covered template sequences. 
If the target ratio cannot be reached, the constraint settings are relaxed up to the relaxation limits.}

\item{fw.primers}{List with optimized primer data frames corresponding to \code{target.temps}. 
Only required for optimizing both strand directions and only 
in the second optimization run in order to check for cross dimerization.}

\item{cur.results.loc}{Directory for storing results of the primer design procedure.}

\item{primer.df}{A data frame of evaluated primer candidates that can be optimized directly.}

\item{updateProgress}{Shiny progress callback function.}
}
\value{
A list containing the results of the primer design procedure:
\describe{
\item{\code{opti}:}{A \code{Primers} object representing the set of optimized primers.}
\item{\code{all_results}:}{A list containing the optimal results for each sampled melting temperature range
                           in terms of a \code{Primers} object 
                           in case that the \code{melting_temp_diff} constraint was active. Otherwise,
                           \code{all_results} only has a single entry representing a primer set relating to an undefined melting temperature.}
\item{\code{used_constraints}:}{A \code{DesignSettings} object with the (adjusted) analysis settings.}
\item{\code{filtered}:}{A \code{Primers} object containing the primer candidates that passed the filtering procedure
                        and which gave rise to the final optimal set.}
}
}
\description{
Designs primers for a single direction.
}
\keyword{internal}

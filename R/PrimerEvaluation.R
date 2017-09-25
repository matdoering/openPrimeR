#' Primer Evaluation.
#' 
#' @rdname PrimerEval
#' @name PrimerEval
#'
#' @description
#' \describe{
#' \item{\code{check_constraints}}{Determines whether a set of primers
#' fulfills the constraints on the properties of the primers.}
#' \item{\code{check_restriction_sites}}{Checks a set of primers 
#' for the presence of restriction sites. 
#' To reduce the number of possible restriction sites,
#' only unambiguous restriction sites are taken into account and 
#' only common (typically used) restriction sites are checked if a common
#' restriction site can be found in a sequence.}
#' \item{\code{filter_primers}}{Filters a primer set according to 
#' the specified constraints such that all primers
#' that do not fulfill the constraints are removed from the primer set.}
#' \item{\code{primer_significance}}{Uses Fisher's exact test to determine 
#' the significance of a primer set according to its ratio of fulfilled 
#' constraints.}
#' \item{\code{subset_primer_set}}{Determines subsets of the input primer set
#' that are optimal with regard to the number of covered template sequences.}
#' }
#'
#' @param primer.df A \code{Primers} object containing the primers
#' whose properties are to be checked.
#' @param template.df A \code{Templates} object containing the 
#' template sequences corresponding to \code{primer.df}.
#' @param settings A \code{DesignSettings} object containing the 
#' constraints that are to be considered.
#' @param active.constraints A subset of the constraint identifiers 
#' provided by \code{settings} that are to be checked
#' for fulfillment. By default \code{active.constraints} is \code{NULL} such that
#' all constraints found in \code{settings} are evaluated. Otherwise,
#' only the constraints specified via \code{active.constraints} 
#' that are available in \code{settings} are considered.
#' @param to.compute.constraints Constraints that are to be computed.
#' By default, \code{to.compute.constraints} is set to \code{NULL} such that
#' all \code{active.constraints} are computed. If \code{to.compute.constraints}
#' is a subset of \code{active.constraints}, all constraints specified
#' via \code{active.constraints} are evaluated for fulfillment,
#' but only the constraints in \code{to.compute.constraints} are newly calculated.
#' @param for.shiny Whether the output of the function shall be
#' formatted as HTML. The default setting is \code{FALSE}.
#' @param updateProgress Progress callback function for shiny. The defaut is
#' \code{NULL} meaning that no progress is monitored via the Shiny interface.
#' @return \code{check_constraints} returns a \code{Primers} object 
#' that is augmented with columns providing the results for the evaluated 
#' constraints.
#' The \code{constraints_passed} column indicates whether all 
#' \code{active.constraints} were fulfilled.
#' The \code{EVAL_*} columns indicate the fulfillment of primer-specific constraints.
#' The \code{T_EVAL_*} columns indicate the fulfillment of template-specific
#' (e.g. coverage-based) constraints.
#' For the coverage computations, columns prefixed by \code{Basic_},
#' indicate the results from string matching, while all other results
#' (e.g. \code{primer_coverage}) indicate the expected coverage
#' after applying the coverage constraints specified in \code{settings}.
#' Columns prefixed by \code{Off_} indicate off-target binding results.
#' @param adapter.action The action to be performed when adapter sequences
#' are found. Either "warn" to issue a warning about adapter sequences or
#' "rm" to remove identified adapter sequences. Currently, only
#' the default setting ("warn") is supported.
#' @param selected Names of restriction sites that are to be checked.
#' By default \code{selected} is \code{NULL} in which case all REBASE 
#' restriction sites are taken into account.
#' @param only.confident.calls Whether only confident calls
#' of restriction sites are returned.
#' All restriction site call is considered \emph{confident} if the restriction site
#' is located in a region that does not match the template sequences.
#' Note that this classification requires that the provided primers
#' are somehow complementary to the provided templates.
#' In contrast, non-confident restriction site calls are 
#' based solely on the primer sequences and do not take the templates
#' into account, resulting in more false positive calls of restriction sites.
#' @param set.name An identifier for the input primers. If \code{NULL},
#' the run identifier is used.
#' @param k The spacing between generated primer subset sizes. By default,
#' \code{k} is set to 1 such that all primer subsets are constructed.
#' @param groups The identifiers of template groups according to which 
#' coverage should be determined. By default, \code{groups} is set to 
#' \code{NULL} such that all all covered templates are considered.
#' @param identifier An identifier for storing the primer set. By default,
#' \code{identifier} is set to  \code{NULL}.
#' @param cur.results.loc Directory for storing the results. By default,
#' \code{cur.results.loc} is set to \code{NULL}, which means that
#' no results are stored.
#' 
#' @note Please note that some constraint computations 
#' may require the installation of additional programs; for more information
#' please view the documentation of \code{\link{DesignSettings}}.
NULL


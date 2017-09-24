#' Coverage Statistics
#'
#' @rdname CoverageStats
#' @name CoverageStats
#'
#' @description
#' \describe{
#' \item{\code{get_cvg_ratio}}{Determines the ratio of template sequences 
#' that are covered by the evaluated input primers. The ratio
#' is in the interval [0,1] where 0 indicates 0\% coverage (no templates
#' covered) and 1 indicates 100\% coverage (all templates covered).}
#' \item{\code{get_cvg_stats}}{Retrieve statistics on covered templates, 
#' either for a single or multiple primer sets.}
#' \item{\code{get_cvg_stats_primer}}{Creates a table summarizing the 
#' coverage events of individual primers.}
#' }
#' 
#' @param primer.df A \code{Primers} object containing the primers.
#' @param template.df A \code{Templates} object containing
#' the template sequences corresponding to \code{primer.df}. 
#' @param allowed.mismatches The number of allowed mismatches
#' for determining the coverage of the templates. By default,
#' all annotated coverage events are considered.
#' @param mode.directionality If \code{mode.directionality} is provided,
#' the coverage of templates is computed for a specific direction of primers.
#' Either "fw" (forward coverage only), "rev" (reverse coverage only), or "both" for both directions. By default, \code{mode.directionality} is \code{NULL}
#' such that the directionality of the primers is determined automatically.
#' @param as.char Whether the coverage ratio should
#' be outputted as a percentage-formatted character vector. By default,
#' \code{as.char} is set to \code{FALSE} such that a numeric is returned.
#' @param primers To retrieve coverage statistics for a single primer set, 
#' please provide an object of class \code{Primers} containing primers with evaluated coverage.
#' To retrieve coverage statistics for multiple primer sets, pelase provide
#' a list with evaluated \code{Primers} objects.
#' @param templates If \code{primers} is an object of class \code{Primers},
#' please provide an object of class \code{Templates} containing the
#' template sequences targeted by \code{primers}. If \code{primers} is a list,
#' \code{templates} should be a list of \code{Template} objects.
#' @param for.viewing Whether the table should be formatted
#' to be human-readable. By default, \code{for.viewing} is \code{FALSE}.
#' @param total.percentages Whether group coverage percentages
#' should be computed in relation to the total number of template sequences
#' or in relation to the number of templates belonging to a specific group.
#' By default, \code{total.percentages} is \code{FALSE} suc that the
#' percentages are group-specific.
#' @param cvg.definition If \code{cvg.definition} is set to
#' "constrained", the statistics for the expected
#' coverage (after applying the coverage constraints) are retrieved.
#' If \code{cvg.definition} is set to "basic", the coverage is determined 
#' solely by string matching (i.e. without applying the coverage constraints).
#' By default, \code{cvg.definition} is set to "constrained".
NULL

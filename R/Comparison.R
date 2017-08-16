#' Evaluated Primer Data for Comparison.
#'
#' Evaluated primer sets targeting the functional human IGH immunoglobulin
#' genes. The sets were generated using the default evaluation settings
#' of openPrimeR. The primer sets were gathered from
#' IMGT and the literature.
#'
#' @docType data
#' @name Comparison
#' @usage data(Comparison)
#'
#' @format \code{primer.data} and \code{template.data} are
#' lists of \code{Primers} and \code{Templates} objects, respectively.
#'
#' @aliases primer.data template.data
#' @keywords datasets
#'
#' @references IMGT®, the international ImMunoGeneTics information system® http://www.imgt.org (founder and director: Marie-Paule Lefranc, Montpellier, France).
#' @examples
#' # Load 'primer.data' and 'template.data'
#' data(Comparison)
#' # Explore the first entry of the primer and template data:
#' primer.data[[1]]
#' template.data[[1]]
#' # Summarize the primer properties:
#' get_comparison_table(template.data, primer.data)
NULL

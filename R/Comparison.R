#' @rdname Data
#' @name Comparison
#' @usage data(Comparison)
#'
#' @aliases primer.data template.data
#'
#' @references IMGT®, the international ImMunoGeneTics information system® http://www.imgt.org (founder and director: Marie-Paule Lefranc, Montpellier, France).
#' @examples
#'
#' # Load the comparison data
#' data(Comparison)
#' # Explore the first entry of the primer and template data:
#' primer.data[[1]]
#' template.data[[1]]
#' # Summarize the primer properties:
#' get_comparison_table(template.data, primer.data)
NULL

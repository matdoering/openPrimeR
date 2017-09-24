#' Template Functionalities
#'
#' @rdname Templates
#' @name Templates
#'
#' @description
#' \describe{
#' \item{\code{adjust_binding_regions}}{Adjusts the existing annotation 
#' of binding regions by specifying
#' a new binding interval relative to the existing binding region.}
#' \item{\code{assign_binding_regions}}{Assigns the primer target binding 
#' regions to a set of template sequences.}
#' \item{\code{update_template_cvg}}{Annotates the template coverage.}
#' \item{\code{select_regions_by_conservation}}{Computes Shannon entropy
#' for the defined binding regions and determines the most conserved regions.}
#' }
#'
#' @param region.fw Interval of new binding regions relative to the forward binding region defined in \code{template.df}.
#' @param region.rev Interval of new binding regions relative to the reverse binding region defined in \code{template.df}
#' @param template.df An object of class \code{Templates}.
#' @param fw Binding regions for forward primers. Either a numeric interval indicating a uniform
#' binding range relative to the template 5' end or a path to a FASTA file providing
#' binding sequences for every template. If \code{fw} is missing, only
#' \code{rev} is considered.
#' @param rev Binding regions for reverse primers. Either a numeric interval indicating a uniform
#' binding range relative to the template 3' end or the path to a FASTA file providing
#' binding sequences for every template. If \code{rev} is missing,
#' only \code{fw} is considered.
#' @param optimize.region If \code{TRUE}, the binding regions
#' specified via \code{fw} and \code{rev} are 
#' adjusted such that binding regions that may form secondary structures are 
#' avoided. This feature requires ViennaRNA (see notes). If \code{FALSE}
#' (the default), the input binding regions are not modified.
#' @param primer.length A numeric scalar providing the probe length that is used for
#' adjusting the primer binding regions when \code{optimize.region} is \code{TRUE}.
#' @param gap.char The character in the input file representing gaps. 
#' @param primer.df An object of class \code{Primers} containing
#' primers with annotated coverage that are to be used to update 
#' the template coverage in \code{template.df}.
#' @param mode.directionality The directionality of primers/templates.
#' @param win.len The extent of the desired primer binding region.
#' This should be smaller than the \code{allowed.region}. The default is 30.
#' @param by.group Shall the determination of binding regions be stratified
#' according to the groups defined in \code{template.df}. By default,
#' this is set to \code{TRUE}.
NULL



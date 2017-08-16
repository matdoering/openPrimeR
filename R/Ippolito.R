#' Evaluated Primer Data from Ippolito et al.
#'
#' Primer and template data for IGHV from Ippolito et al.
#'
#' @docType data
#'
#' @usage data(Ippolito)
#' @name Ippolito
#' @format \code{primer.df} provides a \code{Primers} object containing
#' the evaluated set of primers from Tiller et al.
#' \code{template.df} provides
#' a \code{Templates} object containing functional, human IGHV
#' templates for, and \code{settings} provides a 
#' \code{DesignSettings} object providing the used analysis settings.
#'
#' @keywords datasets
#'
#' @aliases template.df primer.df settings
#'
#' @references Ippolito GC, Hoi KH, Reddy ST, Carroll SM, Ge X, Rogosch T, 
#' Zemlin M, Shultz LD, Ellington AD, VanDenBerg CL, Georgiou G. 2012. 
#' Antibody Repertoires in Humanized NOD-scid-IL2R gamma null Mice and Human B
#' Cells Reveals Human-Like Diversification and Tolerance Checkpoints 
#' in the Mouse. PLoS One 7:e35497.
#' @examples
#' data(Ippolito)
#' # Explore the data:
#' primer.df
#' template.df
#' constraints(settings)
NULL

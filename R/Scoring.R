#' Scoring Functions.
#'
#' @rdname Scoring
#' @name Scoring
#' 
#' @description
#' \describe{
#' \item{\code{score_degen}}{Determines the degeneration score of a sequence.}
#' \item{\code{score_conservation}}{Determines the sequence conservation
#' scores of a set of templates using Shannon entropy.}
#' \item{\code{score_primers}}{Computes scores for a set of primers
#' based on the deviations of the primers from the constraints.}
#' }
#' 
#' @param seq A list of vectors containing individual characters of a nucleotide sequence.
#' @param gap.char The gap character in the sequences.
#' The default is "-".
#' @param template.df A \code{Templates} object providing the set of templates.
#' @param win.len The size of a window for evaluating conservation.
#' The default window size is set to 30.
#' @param by.group Whether the determination of binding regions 
#' should be stratified according to the groups defined in \code{template.df}.
#' The default is \code{TRUE}.
#' @param primer.df A \code{Primers} object containing the primers.
#' @param settings A \code{DesignSettings} object containing the
#' analysis settings.
#' @param active.constraints A character vector of constraint identifiers
#' that are considered for scoring the primers.
#' @param alpha A numeric that is used to determine the trade-off
#' between the impact of the maximal observed deviation and the total
#' deviation. At its default \code{alpha} is set to 0.5 such that
#' the maximal deviation and the total deviation have an equal weight
#' when computing the penalties.
NULL

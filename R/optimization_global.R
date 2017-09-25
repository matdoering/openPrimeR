
penalize_primer <- function(deviations, alpha = 0.5) {
    # alpha: tuning parameter for penalizing maximum deviation
    if (!is.numeric(alpha) || alpha <0 || alpha > 1) {
        stop("Please provide a numeric in the interval [0,1] for 'alpha'.")
    }
    abs.dev <- abs(deviations)
    max.idx <- which.max(abs.dev)
    # minimal penalty is 0. 
    score <- (alpha * abs.dev[max.idx]) + (1 - alpha) * sum(abs.dev)
    return(score)
}
#' @rdname Scoring
#' @name Scoring
#' @details
#' \code{score_primers} determines the penalty of a primer in the following way.
#' Let \code{d} be a vector indicating the absolute deviations from 
#' individual constraints and let \code{p} be the scalar penalty that
#' is assigned to a primer. We define
#' \deqn{p = \alpha \cdot \max_i d_i + \sum_i (1 - \alpha) \cdot d_i}
#' such that for large values of \code{alpha} the maximal deviation 
#' dominates giving rise to a local penalty (reflecting the largest
#' absolute deviation) and for small \code{alpha} the total deviation
#' dominates giving rise to a global penalty 
#' (reflecting the sum of constraint deviations). 
#' When \code{alpha} is 1 only the most extreme absolute deviation is
#' considered and when \code{alpha} is 0 the sum of all absolute 
#' deviations is computed.
#'
#' @return \code{score_primers} returns a data frame containing
#' scores for individual primers.
#' @export
#' @examples
#'
#' # Score the primers
#' data(Ippolito)
#' primer.scores <- score_primers(primer.df, settings)
score_primers <- function(primer.df, settings, active.constraints = names(constraints(settings)), alpha = 0.5) {
    if (!is(primer.df, "Primers")) {
        stop("Please input a 'Primers' object as 'primer.df'.")
    }
    my_penalty <- function(deviation) {
        # local function where 'alpha' doesn't need to be supplied for ddply
        penalize_primer(deviation, alpha = alpha)
    }
    constraint.settings <- constraints(settings)[active.constraints]
    deviations <- get_constraint_deviation_data(primer.df, constraint.settings)
    penalty.df <- ddply(deviations, c("ID"), 
                              here(summarize), Penalty = my_penalty(substitute(Deviation)), Deviation = sum(abs(substitute(Deviation))))
    return(penalty.df)
}


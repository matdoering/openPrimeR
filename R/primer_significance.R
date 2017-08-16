#' Creation of Fulfilled/Failed Constraint Counts.
#'
#' Creates counts of fullfilled/failed constraints.
#'
#' @param primer.df An evaluated \code{Primers} object.
#' @param eval.cols Evaluation columns in \code{primer.df} to consider.
#' By default (\code{NULL}) all evaluation columns are considered.
#' @return A data frame with the number of fulfilled/failed constraints
#' for \code{primer.df}.
#' @keywords internal
create_fulfilled_counts <- function(primer.df, eval.cols = NULL) {
    # compute a data frame of passed/failed constraints
    primer.df <- asS3(primer.df)
    if (length(eval.cols) != 0) {
        m <- match(eval.cols, colnames(primer.df))
        na.idx <- which(is.na(m))
        if (length(na.idx) != 0) {
            msg <- paste("EVAL columns unknown: ", 
                paste(eval.cols, collapse = ","), sep = "")
            warning(msg) 
            m <- m[-na.idx]
        }
        primer.df <- cbind("Run" = primer.df$Run, primer.df[,m, drop = FALSE])
    }
    constraint.idx <- grep("^EVAL_", colnames(primer.df))
    if (length(constraint.idx) == 0) {
        return(NULL)
    }
    constraint.names <- colnames(primer.df)[constraint.idx]
    cols <- c("Run", constraint.names)
    count.df <- plyr::ddply(primer.df[,cols], "Run", 
                    plyr::catcolwise(sum))
    # enrich with failures:
    failure.df <- plyr::ddply(primer.df[,cols], "Run", 
                    plyr::catcolwise(function(x) length(x) - sum(x)))
    failure.names <- paste(constraint.names,"_failure", sep = "")
    colnames(failure.df)[colnames(failure.df) %in% constraint.names] <- 
        failure.names
    count.df <- cbind(count.df, failure.df)
    return(count.df)
}
#' Significance of a Primer Set.
#'
#' Uses Fisher's exact test to determine the significance
#' of a primer set according to its ratio of fulfilled 
#' constraints on the primer properties.
#'
#' The significance is computed by comparing
#' the total count of fulfilled and failed constraints
#' with the corresponding counts of primer sets from the literature.
#' Significant p-values indicate primer sets whose rate of constraint 
#' fulfillment is higher compared to the reference sets.
#'
#' @param primer.df An object of class \code{Primers} for which the
#' significance of physicochemical properties shall be determined.
#' @param set.name An identifier for the input primers. If \code{NULL},
#' the run identifier is used.
#' @param active.constraints Identifiers of the constraints contained in 
#' \code{primer.df} to consider.
#' By default (\code{NULL}) all constraints available in \code{primer.df}
#' are considered for determining the significance.
#' @return The p-value of the primer set according to Fisher's exact test.
#' The returned value has the following attributes: 
#' \describe{
#' \item{\code{test}}{The results of the significance test}
#' \item{\code{tab}}{The confusion matrix for Fisher's exact test}
#' \item{\code{constraints}}{The names of the considered constraints}
#' }
#' @family primer functions
#' @export
#' @examples
#' data(Ippolito)
#' p.data <- primer_significance(primer.df, "Ippolito")
#' attr(p.data,"tab") # the confusion matrix
#' attr(p.data, "test") # results from Fisher's test
#' attr(p.data, "constraints") # considered constraints for the test
primer_significance <- function(primer.df, set.name = NULL, active.constraints = NULL) {
    if (length(primer.df) == 0 || nrow(primer.df) == 0) {
        stop("Cannot compute significance for empty sets.")
    }
    if (is.null(set.name)) {
        set.name <- primer.df$Run[1]
    } 
    if (length(set.name) > 1) {
        stop("Set name should have length 1.")
    }
    # evaluate counts of primer.df
    eval.cols <- NULL
    if (length(active.constraints) != 0) {
        eval.cols <- paste0("EVAL_", active.constraints)
        # select only available constraints
        eval.cols <- eval.cols[eval.cols %in% colnames(primer.df)]
    }
    my.data <- create_fulfilled_counts(primer.df, eval.cols)
    if (is.null(my.data)) { # no evaluated constraints
        warning("Cannot compute a p-value since no constraints have been previously evaluated for the input primers.")
        return(NULL)
    }
    # always exclude specificity and coverage as these are problem-specific
    excl.cons <- c("EVAL_primer_specificity", "EVAL_primer_coverage")
    excl.cons <- c(excl.cons, paste(excl.cons, "_failure", sep = ""))
    constraint.names <- colnames(my.data)[grep("^EVAL_", colnames(my.data))]
    constraint.names <- setdiff(constraint.names, excl.cons)
    if (length(constraint.names) == 0) {
        warning("No computed constraint of the provided primer set is suitable for determining set significance.")
        # no constraints to be tested
        return(NULL)
    }
    ok.names <- constraint.names[!grepl("_failure", constraint.names)]
    failure.names <- constraint.names[grepl("_failure", constraint.names)]
    # now, summarize across all constraints for a fisher's exact test
    #   -> Passed constraint count | Failed constraint count
    my.ok <- unlist(apply(my.data[, ok.names, drop = FALSE], 1, function(x) sum(x, na.rm = TRUE)))
    my.fail <- unlist(apply(my.data[, failure.names, drop = FALSE], 1, function(x) sum(x, na.rm = TRUE)))
    # reference sysdata is available within the package automatically: REF.pass.counts, REF.fail.counts
    m.ok <- match(constraint.names, names(REF.pass.counts))
    m.fail <- match(constraint.names, names(REF.fail.counts))
    ref.ok <- unlist(REF.pass.counts[, m.ok[!is.na(m.ok)]])
    ref.fail <- unlist(REF.fail.counts[, m.fail[!is.na(m.fail)]])
    # test primers:
    tab <- matrix(c(my.ok, my.fail, ref.ok, ref.fail), nrow = 2, 
                  ncol = 2,byrow = TRUE)
    rownames(tab) <- c(set.name, "Reference")
    colnames(tab) <- c("Fulfilled", "Failed")
    # alternative: greater -> check only if odds ratio increases using the tested primers
    # test primers should fulfill more constraints and have less failures
    test.result <- stats::fisher.test(tab, alternative = "greater")
    out.constraint.names <- constraint.names[!grepl("_failure", constraint.names)]
    out.constraint.names <- gsub("EVAL_", "", out.constraint.names)
    result <- test.result$p.value
    attr(result, "test") <- test.result
    attr(result, "table") <- tab
    attr(result, "constraints") <- out.constraint.names
    return(result)
}

 

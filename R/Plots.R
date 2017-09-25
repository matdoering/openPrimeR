#' Plotting Functions.
#'
#' @rdname Plots
#' @name Plots
#' @description
#' \describe{
#' \item{\code{plot_cvg_vs_set_size}}{Plots the coverage ratios of the input primer sets 
#' against the size of the sets.}
#' \item{\code{plot_penalty_vs_set_size}}{Plots the penalties of the input primer sets 
#' against the number of primers contained in each set.
#' The penalties are computed using \code{\link{score_primers}}
#' where more information is provided on how to set \code{alpha}.}
#' \item{\code{plot_primer_subsets}}{Visualizes the coverage of optimized primer subsets.}
#' \item{\code{plot_primer}}{Visualizes the binding positions of every primer relative to
#' the target binding region in the corresponding template sequences.}
#' \item{\code{plot_template_cvg}}{Creates a bar plot visualizing the covered templates.}
#' \item{\code{plot_primer_cvg}}{Shows which groups of templates are covered by individual primers.}
#' \item{\code{plot_constraint}}{Shows the distribution of the primer properties.
#' The current constraint settings are indicated with dashed lines.}
#' \item{\code{plot_constraint_fulfillment}}{Visualizes which primers pass the constraints 
#' and which primers break the constraints}
#' \item{\code{plot_cvg_constraints}}{Plots the distribution of the coverage constraint values.}
#' \item{\code{plot_constraint_deviation}}{Plots the deviation of primer properties from the target ranges.}
#' \item{\code{plot_primer_binding_regions}}{Visualizes the number of binding events of the primers
#' with respect to the allowed binding regions in the templates.}
#' \item{\code{plot_conservation}}{Plots the template sequence conservation (range [0,1]) according to
#' the Shannon entropy of the sequences.}
#' }
#'
#' @param primer.df An object of class \code{Primers} containing
#' primers with evaluated primer coverage.
#' @param template.df An object of class \code{Templates} containing the
#' template sequences.
#' @param primer.data List with objects of class \code{Primers} containing
#' the primer sets that are to be compared.
#' @param primer.subsets A list with optimal primer subsets, each of class \code{Primers}. The \code{k}-th list entry should correspond to an object of class \code{Primers}
#' representing the primer subset of size \code{k} whose coverage ratio
#' is the largest among all possible subsets of size \code{k}.
#' @param template.data List with objects of class \code{Templates} containing
#' the templates corresponding to \code{primer.data}.
#' @param primers Either a single \code{Primers} object with evaluated primer coverage
#' or a list containing such \code{Primers} objects.
#' @param templates If \code{primers} is a \code{Primers} object, \code{templates} should be a \code{Templates} object.
#' If \code{primers} is a list, then \code{templates} should be a list of \code{Templates} objects.
#' @param per.mismatch A logical specifying whether the visualization should be stratified
#' according to the allowed number of mismatches. By default,
#' \code{per.mismatch} is set to \code{FALSE} such that the overall coverage
#' is plotted.
#' @param settings An object of class \code{DesignSettings} containing
#' the constraints to be considered.
#' @param active.constraints A character vector containing the identifiers
#' to be considered for plotting. By default, \code{active.constraints} is
#' \code{NULL} such that all computed constraints found in \code{settings} are plotted.
#' @param show.labels Whether the identifiers of the primer sets
#' should be annotated in the plot. The default is \code{TRUE}.
#' @param highlight.set A character vector providing the identifiers
#' of primer sets to highlight. By default, \code{highlight.set} is \code{NULL}
#' such that no highlighting takes place.
#' @param alpha A numeric in the range [0,1] defining the trade-off between
#' the maximal deviation of a constraint (large \code{alpha}) and
#' all constraint deviations (large \code{alpha}).
#' By default, \code{alpha} is set to 0 such that the absolute
#' deviation across all constraints is considered.
#' @param required.cvg The required coverage ratio.
#' The default is 100\%; this value is plotted as a horizontal line.
#' @param identifier Identifiers of primers that are to be considered.
#' If \code{identifier} is set to \code{NULL} (the default), all primers are considered.
#' @param relation Whether binding positions are computed relative to forward ("fw")
#' or reverse ("rev") binding regions. The default is "fw".
#' @param region.names Character vector of length 2 providing the names
#' of the binding and amplification region.
#' @param ... Optional arguments \code{groups} (a character vector of groups to be plotted when \code{primers} is a single primer set), 
#' \code{highlight.set} (the identifier of a primer set to be highlighted when \code{primers} is a list), 
#' \code{ncol} (a numeric indicating the number of facet columns if \code{primers} is a list),
#' \code{deviation.per.primer} (a boolean indicating whether constraint deviations
#' should be plotted per primer rather than per constraint if \code{primers} is a list)
#' @param direction The directionality of primers to be plotted. This can either
#' be "both" to plot primers of both directions (the default), "fw" to plot
#' only forward primers, or "rev" to plot only reverse primers.
#' @param plot.p.vals An optional logical argument indicating whether
#' p-values computed via \code{\link{primer_significance}} should be annotated in the plot. The default is \code{FALSE}.
#' @param group Optional identifiers of template groups for which binding events should
#' be determined. By default, \code{group} is set to \code{NULL} such that
#' all templates are considered.
#' @param entropy.df A data frame with entropies.
#' Each row gives the entropies of a group of related template
#' sequences for all columns of the alignment.
#' @param alignments A list with \code{DNABin} alignment objects
#' corresponding to the groups (rows) in the alignment.
#' @param gap.char The gap char in the alignments. By default,
#' \code{gap.char} is set to "-".
NULL
 

#' Primer Design Functionalities.
#'
#' @rdname PrimerDesign
#' @name PrimerDesign
#'
#' @description
#' \describe{
#' \item{\code{design_primers}}{Designs a primer set maximizing the number
#' of covered templates using the smallest possible number of primers. 
#' The algorithm tries to ensure that the designed set of primers 
#' achieves a coverage ratio not lower than \code{required.cvg}.
#' To this end, the constraints for designing primers may be relaxed.}
#' \item{\code{get_initial_primers}}{Creates a set of primer candidates 
#' based on the input template sequences. This set of primers can 
#' be used to create custom primer design algorithms.}
#' }
#' 
#' @param template.df A \code{Templates} object containing the template
#' sequences with annotated primer target binding regions.
#' @param mode.directionality The template strand for which primers shall be designed.
#' Primers can be designed either for forward strands ("fw"), 
#' for reverse strands ("rev"), or for both strands ("both"). The default setting
#' is "both".
#' @param settings A \code{DesignSettings} object specifying the constraint settings for designing primers.
#' @param init.algo The algorithm to be used for initializing primers.
#' If \code{init.algo} is "naive", then primers are constructed from substrings of the input template sequences.
#' If \code{init.algo} is "tree", phylogenetic trees are used to form degenerate primers whose degeneracy is bounded by \code{max.degen}.
#' This option requires an installation of MAFFT (see notes). The default \code{init.algo} is "naive".
#' @param opti.algo The algorithm to be used for solving the primer set covering problem. 
#' If \code{opti.algo} is "Greedy" a greedy algorithm is used to solve the 
#' set cover problem. If \code{opti.algo} is "ILP" an integer linear 
#' programming formulation is used. The default \code{opti.algo} is "Greedy".
#' @param required.cvg The desired ratio of of covered template sequences. 
#' If the target coverage ratio cannot be reached, the constraint settings
#' are relaxed according to the the constraint limits in order to reach the target coverage. 
#' The default \code{required.cvg} is set to 1, indicating that 100\% of the templates are to be covered.
#' @param timeout Timeout in seconds. Only applicable when \code{opti.algo} is "ILP".
#' The default is \code{Inf}, which does not limit the runtime.
#' @param max.degen The maximal degeneracy of primer candidates. This setting is particularly
#' relevant when \code{init.algo} is set to "tree". The default setting is \code{16}, which means
#' that at most 4 maximally degenerate positions are allowed per primer.
#' @param conservation Restrict the percentile of considered regions according to their conservation.
#' Only applicable for the tree-based primer initialization. At the 
#' default of 1, all available binding regions are considered.
#' @param sample.name An identifier for the primer design task. The default setting is
#' \code{NULL}, which means that the run identifier provided in \code{template.df} is used.
#' @param cur.results.loc Directory for storing the results of the primer design procedure.
#' The default setting is \code{NULL} such that no output is stored.
#' @param primer.df An optional \code{Primers} object. If an evaluated \code{primer.df} is provided,
#' the primer design procedure only optimizes \code{primer.df} and does not perform
#' the initialization and filtering steps. The default is \code{NULL} such that
#' primers are initialized and filtered from scratch.
#' @param updateProgress Shiny progress callback function. The default is \code{NULL}
#' such that no progress is logged.
#' @param primer.length A scalar numeric providing the 
#' target length of the designed primers. The default length 
#' of generated primers is set to \code{18}.
#' @param primer.estimate Whether the number of required primers shall be estimated. By default (\code{FALSE}), the number of required primers is not estimated.
#' @param sample Character vector providing an identifier for the templates.
#' @param primer.lengths Numeric vector of length 2 providing the 
#' minimal and maximal allowed lengths for generated primers.
#' @param allowed.region.definition A character vector providing the definition
#' of region where primers are to be constructed.
#' If \code{allowed.region.definition} is "within", constructed primers lie within the allowed binding region.
#' If \code{allowed.region.definition} is "any", primers overlap with the allowed binding region.
#' The default is "within".
#' 
#' @note Some constraints can only be computed if additional software is installed,
#' please see the documentation of \code{\link{DesignSettings}} for more information.
#' The usage of \code{init.algo = "tree"} requires an installation of
#' the multiple alignment program MAFFT 
#' (http://mafft.cbrc.jp/alignment/software/).
NULL

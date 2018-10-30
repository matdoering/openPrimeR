#' Settings Functionalities.
#'
#' @rdname Settings
#' @name Settings
#' @description
#' \describe{
#' \item{\code{DesignSettings}}{The \code{DesignSettings} class encapsulates all settings 
#' for designing and evaluating primer sets. 
#' Upon loading an XML file, the \code{DesignSettings} class checks whether
#' the defined constraints can be applied by identifying whether 
#' the requirements for external programs are fulfilled. 
#' If the requirements are not fulfilled, the affected constraints 
#' are removed from the loaded \code{DesignSettings} object
#' and a warning is issued.
#' The loaded constraints are automatically ordered according to
#' the option \code{openPrimeR.constraint_order} such that 
#' the runtime of the \code{\link{design_primers}} and \code{\link{filter_primers}}
#' functions is optimized.}
#' \item{\code{constraints}}{Gets the active constraints of the provided
#' \code{DesignSettings} object.}
#' \item{\code{constraints<-}}{Sets the active constraints of the provided
#' \code{DesignSettings} object.}
#' \item{\code{cvg_constraints}}{Gets the coverage constraints of the provided
#' \code{DesignSettings} object.}
#' \item{\code{cvg_constraints<-}}{Sets the coverage constraints of the provided
#' \code{DesignSettings} object.}
#' \item{\code{conOptions}}{Gets the constraint settings of the provided 
#' \code{DesignSettings} object.}
#' \item{\code{conOptions<-}}{Sets the constraint settings of the provided
#' \code{DesignSettings} object.}
#' \item{\code{constraintLimits}}{Gets the constraint limits that are defined in the provided 
#' \code{DesignSettings} object.}
#' \item{\code{constraintLimits<-}}{Sets the constraint limits of the provided 
#' \code{DesignSettings} object.}
#' \item{\code{PCR}}{Gets the PCR conditions that are defined in the provided
#' \code{DesignSettings} object.}
#' \item{\code{PCR<-}}{Sets the PCR conditions that are defined in the provided
#' \code{DesignSettings} object.}
#' \item{\code{ConstraintSettings}}{The \code{ConstraintSettings} class encapsulates the constraints
#' on the physicochemical properties of primers.}
#' \item{\code{CoverageConstraints}}{The \code{CoverageConstraints} class encapsulates the conditions
#' under which the coverage of primers is evaluated.}
#' \item{\code{PCR_Conditions}}{The \code{PCR_Conditions} class encapsulates the PCR conditions
#' for the computation of primer properties.}
#' \item{\code{ConstraintOptions}}{The \code{ConstraintOptions} class encapsulates the options
#' for constraint computations.}
#' \item{\code{parallel_setup}}{Registers the specified number of cores
#' with the parallel backend.}
#' }
#' @slot Input_Constraints A \code{\link{ConstraintSettings}} object specifying the 
#' desired target value ranges for primer properties. 
#' @slot Input_Constraint_Boundaries A \code{\link{ConstraintSettings}} object specifying
#' the limits for relaxing the constraints during the primer design procedure.
#' This slot may contain the same fields as the \code{Input_Constraints} slot,
#' but the specified desired ranges should be at least as general as
#' those specified in the \code{Input_Constraints} slot.
#' @slot Coverage_Constraints A \code{\link{CoverageConstraints}} object specifying
#' the constraints for computing the primer coverage. 
#' @slot PCR_conditions A \code{\link{PCR_Conditions}} object specifying the PCR-related settings. 
#' @slot constraint_settings A \code{\link{ConstraintSettings}} object providing options 
#' for the computation of individual physicochemical properties.
#' 
#' @slot status Named boolean vector indicating 
#' which of the possible constraints are active (\code{TRUE})
#' and which are not (\code{FALSE}).
#' @slot settings 
#' For \code{ConstraintSettings}, a named list containing the settings for
#' the active constraints.
#' The list may contain the following fields:
#' \describe{
#' \item{primer_coverage:}{The required number of covered template sequences per primer.}
#' \item{primer_specificity:}{The required required specificity of 
#'                            primers in terms of a ratio in the interval [0,1].}
#' \item{primer_length:}{The required lengths of primer sequences.}
#' \item{gc_clamp:}{The desired number of GCs at primer 3' termini.}
#' \item{gc_ratio:}{The desired ratio of GCs in primers 
#'                  in terms of numbers in the interval [0,1].}
#' \item{no_runs:}{The accepted length homopolymer runs in a primer.}
#' \item{no_repeats:}{The accepted length of dinucleotide repeats in a primer.}
#' \item{self_dimerization:}{The lowest acceptable free energy [kcal/mol] for the 
#'                           interaction of a primer with itself. The identification
#'                            of self dimers requires the software \emph{OligoArrayAux} (see notes).}
#' \item{melting_temp_range:}{The desired melting temperature (Celsius) of primers.
#'                           The accurate computation of melting temperatures requires the software \emph{MELTING} (see notes).}
#' \item{melting_temp_diff:}{The maximal allowed difference between the melting temperatures (Celsius)
#'                           of primers contained in the same set. The accurate computation of melting temperatures
#'                           requires the software \emph{MELTING} (see notes).}
#' \item{cross_dimerization:}{The lowest acceptable free energy [kcal/mol] for the
#'                            interaction of a primer with another primer. The identification
#'                            of cross dimers requires the software \emph{OligoArrayAux} (see notes).}
#' \item{secondary_structure:}{The lowest acceptable free energy [kcal/mol] for the
#'                             formation of primer secondary structures. Secondary structures are determined
#'                             using the software \emph{ViennaRNA} (see notes).}
#' }
#' 
#' For \code{PCR_Conditions}, a named list with PCR conditions.
#' The following fields are possible: 
#' \describe{
#' \item{\code{use_taq_polymerase}:}{A logical identifying whether you are performing PCR with a Taq polymerase (\code{TRUE}) or not (\code{FALSE}).}
#' \item{\code{annealing_temp}:}{The annealing temperature in Celsius that is to be used for evaluating the
#' constraints defined in the \code{\link{ConstraintSettings}} object.
#' If the annealing temperature field
#' is not provided, a suitable annealing temperature is automatically computed using a rule of thumb (i.e. subtracting 5 from the melting temperature).}
#' \item{\code{Na_concentration}:}{The molar concentration of monovalent sodium ions.}
#' \item{\code{Mg_concentration}:}{The molar concentration of divalent magnesium ions.}
#' \item{\code{K_concentration}:}{The molar concentration of monovalent potassium ions.}
#' \item{\code{Tris_concentration}:}{The molar concentration of divalent Tris(hydroxymethyl)-aminomethan ions.
#' Note that the Tris ion concentration is about half the buffer concentration.}
#' \item{\code{primer_concentration}:}{The molar concentration of the PCR primers.}
#' \item{\code{template_concentration}:}{The molar concentration of the PCR templates.}
#' }
#' 
#' For \code{CoverageConstraints}, a named list with constraint options. Each 
#' list entry should have an entry \code{min} and/or \code{max}
#' in order to indicate the minimal and maximal allowed values,
#' respectively. 
#' The following identifiers can be used as coverage constraints:
#' \describe{
#' \item{\code{primer_efficiency}:}{The desired efficiencies of primer-template amplification events 
#' in order to be considered as \emph{covered}. \code{primer_efficiency} provides a value in the interval [0,1],
#' which is based on \pkg{DECIPHER}'s thermodynamic model, which considers the impact of 3' terminal mismatches.}
#' \item{\code{annealing_DeltaG}:}{The desired free energies of annealing for putative
#' coverage events between primers and templates. Typically, one would 
#' limit the maximally allowed free energy.}
#' \item{\code{stop_codon}:}{Whether coverage events introducing
#' stop codons into the amplicons should be allowed or discarded. 
#' Here, a value of \code{1} indicates coverage events that induce stop codons.
#' As such, setting both minimum and maximium to zero will disregard
#' coverage events inducing stop codons, while setting the minimum to zero
#' and the maximum to 1 will allow coverage events that induce stop codons.}
#' \item{\code{substitution}:}{Whether coverage events introducing
#' substitutions into the amino acid sequence are considered or discarded.
#' The same encoding as for \code{stop_codon} is used, that is,
#' the value \code{1} indicates coverage events
#' inducing substitutions. Hence, to prevent substitutions,
#' the maximal value of \code{substitution} can be set to zero.}
#' \item{\code{terminal_mismatch_pos}:}{The position relative to 
#' the primer 3' terminal end for which mismatch binding events should be allowed,
#' where the last base in a primer is indicated by position \code{1}.
#' For example, setting the minimal value of \code{terminal_mismatch_pos} 
#' to \code{7} means that only coverage events that do not have a terminal mismatch
#' within the last 6 bases of the primer are allowed.}
#' \item{\code{coverage_model}:}{Use a logistic regression model combining the free energy of annealing and 3' terminal mismatch positions
#' to determine the expected rate of false positive coverage calls. 
#' Using \code{coverage_model}, you can specify the allowed ratio of falsely predicted coverage events.
#' Typically, one would limit the maximal allowed rate of false positives. Note that setting a
#' small false positive rate will reduce the sensitivity of the coverage calls (i.e. true positives will be missed).}
#' }
#' 
#' For \code{ConstraintOptions}, a named list with constraint options.
#' The following fields are permissible:
#' \describe{
#' \item{allowed_mismatches:}{The maximal number of allowed mismatches between
#' a primer and a template sequence. If the number of mismatches of a primer
#' with a template exceeds the specified value, the primer is not considered
#' to cover the corresponding template when the coverage is being computed.}
#' \item{allowed_other_binding_ratio:}{Ratio of allowed binding events
#' outside the target binding ratio. This value should be in the interval
#' [0,1]. If the specified value is greater than zero, all coverage events
#' outside the primer binding region are reported. If, however, the
#' identified ratio of off-target events should exceed the allowed ratio,
#' a warning is issued. If \code{allowed_other_binding_ratio} is set to \code{0},
#' only on-target primer binding events are reported.
#' The setting of \code{allowed_other_binding_ratio} is ignored when designing primers, 
#' which always uses a value of 0.}
#' \item{allowed_region_definition:}{The definition of the target
#' binding regions that is used for evaluating the coverage.
#' In case that \code{allowed_region_definition} is \code{within}, primers have to lie within the allowed binding region.
#' If \code{allowed_region_definition} is \code{any}, primers only have to overlap with the target binding region.}
#' \item{hexamer_coverage:}{If \code{hexamer_coverage} is set to "active", primers whose 3' hexamer (the last 6 bases) is fully complementary to the corresponding
#' template region are automatically considered to cover the template. 
#' If \code{hexamer_coverage} is set to \code{inactive}, 
#' hexamer complementarity does not guarantee template coverage.}
#' }
#'
#' @param x A \code{DesignSettings} object.
#' @param value 
#' An object to be used in one of the setters.
#' For \code{constraints<-} and \code{constraintLimits<-}, a list with constraint settings or boundaries. Each list entry
#' should have a permissible name and consist of at most two
#' values providing the minimal and/or maximal allowed values, which
#' have to be denominated via \code{min} and \code{max}.
#' 
#' For \code{conOptions<-}, a list with constraint options. The permissible
#' fields of the list and their types are documented in the 
#' \code{\link{ConstraintOptions}} class.
#'
#' For \code{cvg_constraints<-}, a list with coverage constraints. Each
#' list entry must have a permissible name and
#' contain a numeric vector with at most two components
#' describing the minimal and/or maximal required values that
#' are to be indicated via \code{min} and \code{max}.
#' The permissible contraint identifiers are documented in the
#' \code{\link{CoverageConstraints}} class.
#' 
#' For \code{PCR<-}, a named list providing PCR conditions 
#' The permissible fields of the list and their types
#' are documented in the \code{\link{PCR_Conditions}} class.
#' @param cores A numeric providing the number of cores to use. The default is \code{NULL}
#' such that half the number of available cores are used.
NULL

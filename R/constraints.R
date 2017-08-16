#' AbstractClass for Constraint Settings.
#'
#' The \code{ConstraintSettings} class encapsulates the constraints
#' on the physicochemical properties of primers.
#'
#' @slot status Named boolean vector indicating 
#' which of the possible constraints are active (\code{TRUE})
#' and which are not (\code{FALSE}).
#' @slot settings Named list containing the 
#' settings of the active constraints
#' @return An \code{AbstractConstraintSettings} object.
#' @keywords internal
AbstractConstraintSettings <- setClass("AbstractConstraintSettings",
	slots = c(status = "vector", settings = "list"),
)

#setMethod("initialize", "AbstractConstraintSettings",
    #function(.Object, ...) {
        #.Object <- callNextMethod(...)
        #.Object
    #}
#)

#' Class for Constraint Settings.
#'
#' The \code{ConstraintSettings} class encapsulates the constraints
#' on the physicochemical properties of primers.
#'
#' You can check whether the constraint settings are fulfilled using
#' \code{\link{check_constraints}} or filter accordingly using
#' \code{\link{filter_primers}}.
#'
#' @section primer_coverage:
#' Computing the primer coverage involves identifying
#' which templates are expected to be amplified (covered) by which primers.
#' The \code{primer_coverage} constraint
#' determines the minimal and maximal number of coverage events per primer
#' that are required. The computation of primer coverage is governed
#' by the coverage constraints postulated via \code{\link{CoverageConstraints}}
#' and the constraint options defined via \code{\link{ConstraintOptions}}. 
#'
#' @section primer_specificity:
#' Primer specificity is automatically determined during the primer coverage
#' computations but the constraint is only checked when the \code{primer_specificity}
#' field is available. The specificity of a primer is defined as its ratio of
#' on-target vs total coverage events (including off-target coverage). Low-specificity
#' primers should be excluded as they may not amplify the target region effectively.
#' 
#' @section primer_length:
#' The length of a primer is defined by its number of bases. Typical
#' primers have lengths between 18 and 22. Longer primers may guarantee higher
#' specificities.
#'
#' @section gc_clamp:
#' The GC clamp refers to the presence of GCs at the 3' end of a primer.
#' For the \code{gc_clamp} constraint, we consider the number of 3' terminal GCs.
#' For example, the primer \emph{actgaaatttcaccg} has a GC clamp of length 3.
#' The presence of a GC clamp is supposed to aid the stability of the
#' polymerase complex. At the same time, long GC clamps should be avoided.
#'
#' @section no_runs:
#' Homopolymer runs (e.g. the primer \emph{aaaaa} has a run of 5 A's)
#' may lead to secondary structure formation and unspecific binding
#' and should therefore be avoided.
#'
#' @section no_repeats:
#' Dinucleotide repeats (e.g. the primer \emph{tatata} has 3 TA repeats)
#' should be avoided for the same reason a long homopolymer runs.
#'
#' @section self_dimerization:
#' Self dimerization refers to a primer that binds to itself rather than
#' to one of the templates. Primers exhibiting self dimers should
#' be avoided as they may prevent the primer from amplifying the templates.
#' Therefore primers with small free energies of dimerization should 
#' be avoided.
#'
#' @section melting_temp_range:
#' The melting temperature is the temperature at which 50% of the primers
#' are in duplex with templates and 50% are still present as random coils.
#' Hence, primers exhibiting high melting temperatures have high affinities
#' to the templates, while primers with small melting temperatures
#' have small affinities. The melting temperatures of the primers
#' determine the annealing temperature of the PCR, which is why
#' the melting temperatures of the primers should not deviate too much (see
#' \code{melting_temp_diff}).
#'
#' @section melting_temp_diff:
#' The differences between the melting temperatures of primers in a
#' set of primers should not deviate too much as the annealing temperaturte
#' of a PCR should be based on the smallest melting temperature
#' of a primer in the set. If there are other primers in the set 
#' exhibiting considerably higher melting temperatures,
#' these primers may bind inspecifically due to the low annealing temperature.
#'
#' @section cross_dimerization:
#' When two different primers bind to each each other rather than
#' to the templates, this is called cross dimerization.
#' Cross dimerization should be prevent at all costs because
#' such primers cannot effectively amplify their target templates.
#' Cross dimerizing primers can be excluding primers
#' exhibiting small free energies of cross dimerization.
#'
#' @section secondary_structure:
#' When a primer exhibits secondary structure, this may prevent it from
#' binding to the templates. To prevent this, 
#' primers with low free energies of secondary structure formation
#' can be excluded.
#'
#' @slot status Named boolean vector indicating 
#' which of the possible constraints are active (\code{TRUE})
#' and which are not (\code{FALSE}).
#' @slot settings Named list containing the 
#' settings of the active constraints:
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
#' @return A \code{ConstraintSettings} object.
#' @note
#' External programs are required for computing the following constraints:
#' \describe{
#' \item{MELTING (http://www.ebi.ac.uk/biomodels/tools/melting/):}{Thermodynamic computations (optional) for determining melting temperatures for the constraints \code{melting_temp_diff} and \code{melting_temp_range}}
#' \item{OligoArrayAux (http://unafold.rna.albany.edu/OligoArrayAux.php):}{Thermodynamic computations used for computing \code{self_dimerization} and \code{cross_dimerization}.
#' Also required for computing \code{primer_coverage} when a constraint based on the free energy of annealing is active.}
#' \item{ViennaRNA (http://www.tbi.univie.ac.at/RNA/):}{Secondary structure predictions used for the constraint \code{secondary_structure}}
#' }
#' @family settings functions
#' @export
#' @keywords Classes
#' @examples
#' # Initializing a new 'ConstraintSettings' object:
#' constraint.settings <- new("ConstraintSettings")
#' # Retrieving the constraint settings from a 'DesignSettings' object:
#' data(Ippolito) # loads a 'DesignSettings' object into 'settings'
#' constraints(settings)
#' # Modifying the constraint settings:
#' constraints(settings)$no_runs["max"] <- 10
#' constraints(settings) <- constraints(settings)[names(constraints(settings)) != "gc_clamp"]
ConstraintSettings <- setClass("ConstraintSettings", contains="AbstractConstraintSettings",
	slots = c( # constructor arguments:
			status = "vector",
            settings = "list"
			),
    # define initial values using prototype
    prototype = prototype(status = {
            possible.constraints <- c("primer_coverage", "primer_length", 
                "primer_specificity", "gc_clamp", "gc_ratio","no_runs",
                       "no_repeats", "self_dimerization", "melting_temp_diff",
                       "melting_temp_range", "cross_dimerization", 
                       "secondary_structure")
            cons <- rep(FALSE, length(possible.constraints))
            names(cons) <- possible.constraints
            cons
        },
        settings = list() # no active constraints initially
   ),
    validity=function(object)
	{
		return(check_constraint_settings_validity(object))
	}
)
setMethod("initialize", "ConstraintSettings",
    function(.Object, ...) {
        .Object <- callNextMethod(...)  # passed to parent constructor
        .Object
})
setMethod("show", "ConstraintSettings", function(object) {
    print(constraints.to.df(object@settings, out.names = "Target range"))
    invisible()
})

#' Check the Validity of the Constraint Settings.
#'
#' Checks whether the status and the active constraints match.
#' Determines whether the constraints are allowed/known.
#'
#' @param object An \code{AbstractConstraintSettings} object.
#' @return Error messages in case of errors, otherwise \code{TRUE}.
#' @keywords internal
check_constraint_settings_validity <- function(object) {
    # store errors here:
    errors <- character(0)
    # create reference objects for checking possible fields of the settings
    if (is(object, "ConstraintSettings")) {
        ref.object <- ConstraintSettings() # containing all possible constraints
    } else if (is(object, "CoverageConstraints")) {
        ref.object <- CoverageConstraints()
    } else if (is(object, "ConstraintOptions")) {
        ref.object <- ConstraintOptions()
    } else if (is(object, "PCR_Conditions")) {
        ref.object <- PCR_Conditions()
    } else {
        msg <- paste0("Unknown class: ", class(object))
        errors <- c(errors, msg)
    }
    # possible constraints should always stay the same:
    if (length(names(object@status)) == 0 || 
    names(ref.object@status) != names(object@status)) {
        m <- match(names(object@status), names(ref.object@status))
        idx <- which(is.na(m))
        error <- paste0("'status' was modified from the default.",
                    "Unknown constraint: ",
                    paste0(names(object@status)[idx], collapse = ","))
        errors <- c(errors, error)
    } 
    # if a constraint is active, it should have a setting:
    active.cons <- names(object@status)[which(object@status)]
    if (any(!active.cons %in% names(object@settings))) {
        error <- "All active constraints should have a setting."
        errors <- c(errors, error)
    }
    # if a constraint has a setting, it should be active
    setting.cons <- names(object@settings)
    if (length(setting.cons) != 0 && any(!setting.cons %in% names(object@status))) {
        error <- "Constraint with setting wasnt' active."    
        errors <- c(errors, error)
    }
    # check the types of the constraint settings:
    if (is(object, "ConstraintSettings") || is(object, "CoverageConstraints")) {
        #print("cvg interval check")
        #print(object@settings)
        # check interval properties for the constraints
        classes <- sapply(object@settings, class)
        if (!all(classes %in% c("numeric", "integer"))) {
            msg <- "Please provide only numeric constraints."
            errors <- c(errors, msg)
        }
        # check for interval property of the constraints:
        interval.check <- check_interval(object@settings)
        #print("checked interval:")
        #print(interval.check)
        if (!is.logical(interval.check)) {
            errors <- c(errors, interval.check)
        }
    } else if (is(object, "ConstraintOptions") || is(object, "PCR_Conditions")) {
        # check mandatory options for the other settings
        m <- match(names(ref.object@settings), names(object@settings))
        if (any(is.na(m))) {
            idx <- which(is.na(m))
            msg <- paste0("There were missing mandatory settings: ", 
                 paste0(names(constraints(ref.object))[idx], collapse = ","))
            errors <- c(errors, msg)
        }
        # check values of options
        if (is(object, "ConstraintOptions")) {
            # character options:
            region.ok <- match.arg(constraints(object)$allowed_region_definition, c("within", "any"))
            hexamer.ok <- match.arg(constraints(object)$hexamer_coverage, c("active", "inactive"))
            # other options should be numeric: 
            sel <- which(!names(constraints(object)) %in% c("allowed_region_definition", "hexamer_coverage"))
            classes <- sapply(constraints(object)[sel], class)
            if (!all(classes %in% c("numeric", "integer"))) {
                msg <- paste0("Incorrect option class: should have been numeric.")
                errors <- c(errors, msg)
            }
        }
        # check values of PCR
        if (is(object, "PCR_Conditions")) {
            if (!is(constraints(object)$use_taq_polymerase, "logical")) {
                msg <- "'use_taq_polymerase' should be logical."
                errors <- c(errors, msg)
            }
            # other PCR conditions should be numeric:
            sel <- which(!names(constraints(object)) %in% c("use_taq_polymerase"))
            classes <- sapply(constraints(object)[sel], class)
            if (!all(classes %in% c("numeric", "integer"))) {
                msg <- paste0("Incorrect class for PCR setting: should have been numeric.")
                errors <- c(errors, msg)
            }
        }
    } 
    if (length(errors) == 0) {
        return(TRUE)
    } else {
        return(errors)
    }
}
#' Class for Coverage Constraints.
#'
#' The \code{CoverageConstraints} class encapsulates the conditions
#' under which the coverage is evaluated.
#'
#' @slot status Named boolean vector indicating 
#' which of the possible options are active (\code{TRUE})
#' and which are not (\code{FALSE}).
#' @slot settings Named list with constraint options. Each 
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
#' @return A \code{CoverageConstraints} object.
#' @note
#' External programs are required for computing the following constraints:
#' \describe{
#' \item{OligoArrayAux (http://unafold.rna.albany.edu/OligoArrayAux.php):}{Thermodynamic computations used for 
#' computing the coverage constraints \code{annealing_DeltaG}, \code{primer_efficiency}, and \code{coverage_model}}
#' }
#' @family settings functions
#' @export
#' @keywords Classes
#' @examples
#' # Initialize a new 'CoverageConstraints' object:
#' cvg.constraints <- new("CoverageConstraints")
#' # Retrieving the coverage constraints from a 'DesignSettings' object:
#' data(Ippolito) # loads a 'DesignSettings' object into 'settings'
#' cvg_constraints(settings)
#' # Modifying the coverage constraints
#' cvg_constraints(settings)$primer_efficiency["min"] <- 0.001
CoverageConstraints <- setClass("CoverageConstraints", contains="AbstractConstraintSettings",
    prototype = prototype(status = {
            possible.constraints <- c("coverage_model", "stop_codon", "substitution", 
                "terminal_mismatch_pos", "primer_efficiency", "annealing_DeltaG", 
                "coverage_model")
            cons <- rep(FALSE, length(possible.constraints))
            names(cons) <- possible.constraints
            cons
        },
        settings = list() # no active constraints initially
   ),
    validity=function(object)
	{
		return(check_constraint_settings_validity(object))
	}
)
setMethod("show", "CoverageConstraints", function(object) {
    print(constraints.to.df(object@settings, out.names = "Setting"))
    invisible()
})
setMethod("initialize", "CoverageConstraints",
    function(.Object, ...) {
        .Object <- callNextMethod(...)  # passed to parent constructor
        .Object
})
#' Class for Constraint Options.
#'
#' The \code{ConstraintOptions} class encapsulates the options
#' for constraint computations.
#'
#' @slot status Named boolean vector indicating 
#' which of the possible options are active (\code{TRUE})
#' and which are not (\code{FALSE}).
#' @slot settings Named list with constraint options.
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
#'}
#'
#' @return A \code{ConstraintOptions} object.
#' @family settings functions
#' @export
#' @keywords Classes
#' @examples
#' # Initialize a new 'ConstraintOptions' object:
#' constraint.options <- new("ConstraintOptions")
#' # Retrieve the constraint options from a 'DesignSettings' object:
#' data(Ippolito) # loads a 'DesignSettings' object into 'settings'
#' conOptions(settings)
#' # Prevent off-target binding:
#' conOptions(settings)$allowed_other_binding_ratio <- 0
ConstraintOptions <- setClass("ConstraintOptions", contains="AbstractConstraintSettings",
    prototype = prototype(status = {
        # set mandatory options to TRUE
        option.status <- c("allowed_mismatches" = TRUE, "allowed_other_binding_ratio" = TRUE,
                          "allowed_region_definition" = TRUE, "hexamer_coverage" = FALSE)
        option.status
    }, 
    # ensure that 'settings' reflects the option status (mandatory should occur in settings):
    settings = list("allowed_mismatches" = 7, "allowed_other_binding_ratio" = 1, 
                    "allowed_region_definition" = "within")
   ),
   validity=function(object) {
		return(check_constraint_settings_validity(object))
	}
)
setMethod("show", "ConstraintOptions", function(object) {
    print(create.options.table(object@settings))
    invisible()
})
setMethod("initialize", "ConstraintOptions",
    function(.Object, ...) {
        .Object <- callNextMethod(...)  # passed to parent constructor
        .Object
})
#' Class for PCR Conditions.
#'
#' The \code{PCR_Conditions} class encapsulates the PCR conditions
#' for the computation of primer properties.
#'
#' @slot status Named boolean vector indicating 
#' which of the possible options are active (\code{TRUE})
#' and which are not (\code{FALSE}).
#' @slot settings Named list with constraint options.
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
#' \item{\code{Tris_concentration}:}{The molar concentration of the Tris(hydroxymethyl)-aminomethan buffer. 
#' Note that this value is the buffer concentration. To determine corresponding Tris ion concentration, the value of the buffer concentration is halved.}
#' \item{\code{primer_concentration}:}{The molar concentration of the PCR primers.}
#' \item{\code{template_concentration}:}{The molar concentration of the PCR templates.}
#' }
#' @return A \code{PCR_Conditions} object.
#' @family settings functions
#' @keywords Classes
#' @export
#' @examples
#' # Initialize a new 'PCR_Conditions' object:
#' PCR.conditions <- new("PCR_Conditions")
#' # Retrieving the PCR conditions from a 'DesignSettings' object:
#' data(Ippolito) # loads a 'DesignSettings' object into 'settings'
#' PCR(settings)
#' # Modifying the PCR conditions:
#' PCR(settings)$use_taq_polymerase <- FALSE
PCR_Conditions <- setClass("PCR_Conditions", contains="AbstractConstraintSettings",
    prototype = prototype(status = {
        # set mandatory options to TRUE
        option.status <- c("use_taq_polymerase" = TRUE, "annealing_temp" = FALSE,
                        "Na_concentration" = TRUE, "Mg_concentration" = TRUE,
                        "K_concentration" = TRUE, "Tris_concentration" = TRUE,
                        "primer_concentration" = TRUE, "template_concentration" = TRUE,
                        "cycles" = FALSE)
        option.status
    }, 
    # ensure that 'settings' reflects the option status (mandatory should occur in settings):
    settings = list("use_taq_polymerase" = TRUE, "Na_concentration" = 0,
                    "Mg_concentration" = 0.0015, "K_concentration" = 0.05, 
                    "Tris_concentration" = 0, "primer_concentration" = 2e-07,
                    "template_concentration" = 1.28e-11)
   ),
   validity=function(object) {
		return(check_constraint_settings_validity(object))
	}
)
setMethod("show", "PCR_Conditions", function(object) {
    print(create.PCR.table(object@settings))
    invisible()
})

setMethod("initialize", "PCR_Conditions",
    function(.Object, ...) {
        .Object <- callNextMethod(...)  # passed to parent constructor
        .Object
})

#' Class for Primer Design Settings.
#'
#' The \code{DesignSettings} class encapsulates all settings 
#' for designing and evaluating primer sets. 
#' Upon loading an XML file, the \code{DesignSettings} class checks whether
#' the defined constraints can be applied by identifying whether 
#' the requirements for external programs are fulfilled. 
#' If the requirements are not fulfilled, the affected constraints 
#' are removed from the loaded \code{DesignSettings} object
#' and a warning is issued.
#' The loaded constraints are automatically ordered according to
#' the option \emph{openPrimeR.constraint_order} such that 
#' the runtime of the \code{\link{design_primers}} and \code{\link{filter_primers}}
#' functions is optimized. 
#'
#' Note that the fields \code{Input_Constraints}, \code{Input_Constraint_Boundaries}, and \code{Coverage_Constraints} should 
#' contain entries with at most two components using the fields \code{min} and/or \code{max}.
#'
#' The \code{Input_Constraint_Boundaries} should always be at least as general as the
#' specified \code{Input_Constraints}.
#'
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
#' @name DesignSettings-class
#' @rdname DesignSettings-class
#' @exportClass DesignSettings
#' @return A \code{DesignSettings} object.
#' 
#' @family settings functions
#' @keywords Settings
#' @seealso \code{\link{read_settings}} for reading settings from XML files,
#' \code{\link{write_settings}} for storing settings as XML files,
#' \code{\link{constraints}} for accessing constraints,
#' \code{\link{constraintLimits}} for accessing constraint boundaries,
#' \code{\link{cvg_constraints}} for accessing coverage constraints,
#' \code{\link{conOptions}} for accessing constraint options,
#' \code{\link{PCR}} for accessing the PCR conditions.
#' @examples
#' # Load a settings object
#' filename <- system.file("extdata", "settings", 
#'                  "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
#' settings <- read_settings(filename)
#' # Modify the constraints
#' constraints(settings)$gc_clamp["min"] <- 0
#' # Modify the constraint limits for designing primers
#' constraintLimits(settings)$gc_clamp["max"] <- 6
#' # Modify the coverage constraints
#' cvg_constraints(settings)$primer_efficiency["min"] <- 0.001
#' # Modify the PCR conditions
#' PCR(settings)$Na_concentration <- 0.0001
#' # Modify the constraint options
#' conOptions(settings)$allowed_mismatches <- 0
DesignSettings <- setClass("DesignSettings",
	slots = c( # constructor arguments:
			Input_Constraints = "ConstraintSettings",
			Input_Constraint_Boundaries = "ConstraintSettings",
            Coverage_Constraints = "CoverageConstraints",
			PCR_conditions = "PCR_Conditions",
			constraint_settings = "ConstraintOptions"
			),
	validity=function(object)
	{
		return(check_settings_validity(object))
	}
)
#' Check of limit correctness.
#'
#' Checks whether a constraint limit is more general than the setting.
#'
#' @param setting A single constraint setting.
#' @param limit A single constraint limit.
#' @return A vector containing \code{TRUE} if the limit is more general
#' than the constraint setting and \code{FALSE} otherwise.
#' @keywords internal
check_limit_value <- function(setting, limit) {
    if (length(limit) == 0) { # no limit specified
        # nothing to check
        return(TRUE)
    }
    if (length(names(setting)) == 0 || !all(names(setting)  %in% c("min", "max"))) {
        stop("Please provide min/max as fields for constraints.")
    }
    if (length(names(limit)) == 0 || !all(names(limit) %in% c("min", "max"))) {
        stop("Please provide min/max as fields for constraint limits.")
    }
    ret <- NA
    if (length(setting) == 2) {
        ret <- c(limit[1] <= setting[1], limit[2] >= setting[2])
    } else if (length(setting) == 1 && names(setting) == "max") {
        ret <- limit["max"] >= setting["max"]
    } else if (length(setting) == 1 && names(setting) == "min") {
        ret <- limit["min"] <= setting["min"]
    } else {
        message("Unkown constraint length:")
        message(setting)
        message(limit)
        ret <- NA
    }
    na.idx <- which(is.na(ret))
    ret[na.idx] <- FALSE # missing value in one of the settings
    return(ret)
}
#' Correction of Constraint Boundaries.
#'
#' Fixes the constraint boundaries if they are more narrow than the current settings.
#'
#' @param constraints A list with constraint settings.
#' @param constraint.limits A list with constraint limits.
#' @param fix.limit Whether the constraint limits should be
#' adjusted. If \code{FALSE}, the constraint settings are adjusted.
#' @return The corrected constraint limits.
#' @keywords internal
fix_constraint_boundaries <- function(constraints, constraint.limits, fix.limit = TRUE) {
    iter.range <- NULL
    if (fix.limit) {
        iter.val <- constraints # the constraints for iteration
        other.val <- constraint.limits # the constraints to be adjusted
    } else {
        iter.val <- constraint.limits
        other.val <- constraints
    }
    for (i in seq_along(iter.val)) {
        con <- names(iter.val)[i]
        if (!con %in% names(other.val)) {
            next
        }
        other.setting <- other.val[[con]]
        setting <- iter.val[[con]]
        #print(con)
        #print(setting)
        #print(other.setting)
        if (fix.limit) {
            test <- check_limit_value(setting, other.setting)
        } else {
            test <- check_limit_value(other.setting, setting)
        }
        if (!all(test)) {
            #warning(paste("Had to fix constraint limit of constraint:",
                #con), ". Constraint limits should be more general than the setting.")
            # limit needs to be adjusted to the current setting
            other.setting[!test] <- setting[!test]
            names(other.setting)[!test] <- names(setting)[!test]
            other.val[[con]] <- other.setting
        } 
    }
    # sync the constraints with the constraintlimits: addition/deletion of constraints
    # remove all constraints that are present in 'other.val' but not in 'iter.val'
    rm.cons <- setdiff(names(other.val), names(iter.val))
    if (length(rm.cons) != 0) {
        other.val <- other.val[-which(names(other.val) %in% rm.cons)]
    }
    # add constraints that are present in 'iter.val' but not in 'other.val'
    add.cons <- setdiff(names(iter.val), names(other.val))
    # set missing entries to have the identical constraint and boundary since we don't know better:
    other.val[add.cons] <- iter.val[add.cons]
    # ensure that constraints are ordered correctly 
    # (the order agrees with the order of the options, since we ordered in constraints()/constraintLimits())
    con.order <- names(iter.val)
    other.val <- other.val[con.order[con.order %in% names(other.val)]]
    return(other.val)
}

setMethod("show", "DesignSettings", function(object) {
    con.tab <- create.constraint.table(constraints(object), 
                                       constraintLimits(object))
    # raw constraint name output: users should know the identifiers that are correct from the output.
    inactive.cons <- setdiff(names(object@Input_Constraints@status), names(constraints(object)))
    con.result <- list("Active" = con.tab, "Inactive" = inactive.cons)
    cvg.con.tab <- create.constraint.table(cvg_constraints(object), NULL)
    inactive.cons <- setdiff(names(object@Coverage_Constraints@status), names(cvg_constraints(object)))
    cvg.con.result <- list("Active" = cvg.con.tab, "Inactive" = inactive.cons)
    PCR.tab <- create.PCR.table(PCR(object))
    inactive.options <- setdiff(names(object@PCR_conditions@status), names(PCR(object)))
    PCR.result <- list("Active" = PCR.tab, "Inactive" = inactive.options)
    opt.tab <- create.options.table(conOptions(object))
    inactive.options <- setdiff(names(object@constraint_settings@status), names(conOptions(object)))
    opt.result <- list("Active" = opt.tab, "Inactive" = inactive.options)
    out <- list("Constraints" = con.result, "Coverage Constraints" = cvg.con.result,
                "PCR conditions" = PCR.result, "Options" = opt.result)
    print(out)
})

setMethod("initialize", "DesignSettings",
    function(.Object, constraint.settings, constraint.limits, cvg.constraints,
             PCR.conditions, constraint.options) {
        # just call the default constructor
        #print("INITIALIZATION:")
        # .Object is the prototype from the class for which slots should be set
        if (is(constraint.settings, "ConstraintSettings")) {
            # ok
        } else if (is(constraint.settings, "list")) {
            obj <- ConstraintSettings()
            constraints(obj) <- constraint.settings
            constraint.settings <- obj
        } else {
            stop("Please provide 'constraint.settings' as a 'ConstraintSettings' object or a list that can be converted to a 'ConstraintSettings' object.")
        }
        .Object@Input_Constraints <- constraint.settings
        if (is(constraint.limits, "ConstraintSettings")) {
            # ok
         } else if (is(constraint.limits, "list")) {
            obj <- ConstraintSettings()
            constraints(obj) <- constraint.limits
            constraint.limits <- obj
         } else {
            stop("Please provide 'constraint.limits' as a 'ConstraintSettings' object or a 'list'.")
        }
        .Object@Input_Constraint_Boundaries <- constraint.limits
        if (is(cvg.constraints, "CoverageConstraints")) {
            # ok
        } else if (is(cvg.constraints, "list")) {  
            obj <- CoverageConstraints()
            constraints(obj) <- cvg.constraints
            cvg.constraints <- obj
        } else {
            stop("Please provide 'cvg.constraints' as a 'CoverageConstraints' object or a 'list'.")
        }
        .Object@Coverage_Constraints <- cvg.constraints
        if (is(PCR.conditions, "PCR_Conditions")) {
            # ok
        } else if (is(PCR.conditions, "list")) {
            obj <- PCR_Conditions()
            constraints(obj) <- PCR.conditions
            PCR.conditions <- obj
        } else {
            stop("Please provide 'PCR.conditions' as a 'PCR_Conditions' object or a 'list'.")
        }
        .Object@PCR_conditions <- PCR.conditions
        if (is(constraint.options, "ConstraintOptions")) {
            # ok
        } else if (is(constraint.options, "list")) {
            obj <- ConstraintOptions()
            constraints(obj) <- constraint.options
            constraint.options <- obj
        } else {
            stop("Please provide 'constraint.options' as a 'ConstraintOptions' object.")
        }
        .Object@constraint_settings <- constraint.options
        if (is(.Object, "try-error")) {
            stop("Could not initialize DesignSettings: did you supply all args?\n", attr(.Object, "condition"))
        }
        # specify the order of constraint evaluation for the design procedure 
        # users can modify the order by changing the options() setting
        con.order <- getOption("openPrimeR.constraint_order")
        # reduce possible constraints according to available tools
        sel.constraints <- con_select(names(constraints(.Object@Input_Constraints)))
        rm.cons <- setdiff(names(constraints(.Object@Input_Constraints)), sel.constraints)
        if (length(rm.cons) != 0) {
            message("Due to missing external tools the following constraints were ignored:",
                    paste(rm.cons, collapse = ","))
        }
        sel.o <- con.order[con.order %in% sel.constraints]
        constraints(.Object@Input_Constraints) <- constraints(.Object@Input_Constraints)[sel.o]
        sel.constraints <- con_select(names(constraints(.Object@Input_Constraint_Boundaries)))
        sel.o <- con.order[con.order %in% sel.constraints]
        constraint.boundaries <- fix_constraint_boundaries(constraints(.Object@Input_Constraints), constraints(.Object@Input_Constraint_Boundaries)[sel.o])
        constraints(.Object@Input_Constraint_Boundaries) <- constraint.boundaries
        cvg.sel <- select.constraints(names(constraints(.Object@Coverage_Constraints))) # select.constraints here since we're not storing cvg constraints in option
        rm.cvg <- setdiff(names(constraints(.Object@Coverage_Constraints)), cvg.sel)
        if (length(rm.cvg) != 0) {
            message("Due to missing external tools the following coverage constraints were ignored:",
                    paste(rm.cvg, collapse = ","))
        }
        constraints(.Object@Coverage_Constraints) <- constraints(.Object@Coverage_Constraints)[cvg.sel]
        validObject(.Object)
        .Object
    }
) 

#' Check Setting Names.
#'
#' Checks whether the specified settings hvae the correct names.
#'
#' @param known.options Allowed setting names.
#' @param input.options Input setting names
#'
#' @return Mapping of \code{input.options} to \code{known.options} or \code{NULL} if invalid.
#' @keywords internal
check_names <- function(known.options, input.options) {
    # checks whether the names of 
    errors <- NULL
    m <- match(input.options, known.options)
    if (any(is.na(m))) {
        # unknown option
        errors <- c(errors, paste("Unknown/missing field in object:", paste(input.options[is.na(m)], collapse = ",")))
    }
    if (length(errors) != 0) {
        return(errors)
    } else {
        return(m)
    }
}

#' Check Setting Validity.
#'
#' Checks whether the input settings are valid or not.
#'
#' @param known.options Vector with names and classes of allowed options.
#' @param options Active options to be checked.
#' @param mandatory.options Fields that have to be present.
#' @return \code{TRUE} if the setting is valid, \code{FALSE} otherwise.
#' @keywords internal
check_setting <- function(known.options, options, mandatory.options = NULL) {
    # checks whether a setting in DesignSettings is valid or not.
    if (length(options) == 0) {
        return(TRUE)
    }
    # check for correct field names:
    errors <- NULL
    if (!is.null(known.options)) {
        m <- check_names(names(known.options), names(options))
        if (is.character(m)) { # error
            return(m)
        }
        c <- sapply(seq_along(options), function(x) !(class(options[[x]]) %in% known.options[m][[x]]) & !is.null(options[[x]]))
        if (any(c)) {
            idx <- which(c)
            msgs <- sapply(idx, function(x) paste("Supplied class incorrect: ",  names(options[which(c)]), 
                    ".Class was: ", sapply(options, class)[c], ", but should have been: ", 
                    paste(known.options[m][[which(c)]], collapse = " or "), ".", sep = ""))
            return(msgs)
        }
    }
    # check for mandatory field names:
    if (!is.null(mandatory.options)) {
        m <- check_names(names(options), names(mandatory.options))
        if (is.character(m)) {
            # error
            return(m)
        }
        if (is(options, "data.frame")) {
            c <- sapply(seq_along(mandatory.options), function(x) !(class(asS3(options)[[m[x]]]) %in% mandatory.options[[x]]))
        } else {
            c <- sapply(seq_along(mandatory.options), function(x) !(class(options[m][[x]]) %in% mandatory.options[[x]]))
        }
        if (any(c)) {
            idx <- which(c)
            msgs <- sapply(idx, function(x) paste("Supplied class was incorrect for '",  names(mandatory.options)[x],
                    "'. Class was '", sapply(options, class)[m[x]], "', but should have been '", 
                    paste(mandatory.options[[x]], collapse = " or "), "'.", sep = ""))
            return(msgs)
        }
    }
    return(TRUE)
}

#' Check Constraint Intervals
#'
#' Checks the validity of constraint intervals.
#'
#' @param constraints A list with constraint settings.
#' @return \code{TRUE}, if all constraints specificy valid intervals, 
#' \code{FALSE} otherwise.
#' @keywords internal
check_interval <- function(constraints) {
    check <- rep(FALSE, length(constraints))
    for (i in seq_along(constraints)) {
        con <- constraints[[i]]
        #print("check_interval:")
        #print(names(constraints)[i])
        #print(con)
        if (!all(is.numeric(con) | is.integer(con))) {
            return("Please supply only integer/numeric values.")
        }
        if (length(names(con)) == 0) {
            check[i] <- FALSE
        } else if (length(con) == 1) {
            # only min or max condition
            if (names(con) %in% c("min", "max")) {
                check[i] <- TRUE
            }
        } else if (length(con) == 2) {
            # min & max conditions
            if (names(con) %in% c("min", "max") &&
                con["min"] <= con["max"]) {
                check[i] <- TRUE
            }
        } else {
            # other lengths are not allowed
            check[i] <- FALSE
        }
    }
    if (all(check)) {
        return(TRUE)
    } else {
        details <- names(constraints)[!check]
        msg <- (paste("Constraint interval (min/max) was not specified properly ",
                    "for the following constraints:", 
                    paste(details, collapse = ","), 
                    ".\nDid you supply a named vector containing the components 'min' and/or 'max' such that 'min' <= 'max' if both entries were indeed supplied?", sep = ""))
        return(msg)
    }
}
#' Validity Check for Limits.
#'
#' Checks whether the constraint limits are at least as
#' general as the constraint settings. This ensures that
#' the relaxation works in the proper direction.
#'
#' @param constraint.settings A list with the constraint settings.
#' @param constraint.limits A list with the constraint relaxation limits.
#' @return \code{TRUE} if the limits are at least as wide as the constraints, 
#' \code{FALSE} otherwise.
#' @keywords internal
check_limits <- function(constraint.settings, constraint.limits) {
    check <- sapply(seq_along(constraint.settings), function(x) {
        id <- names(constraint.settings)[x]
        setting <- constraint.settings[[id]]
        limit <- constraint.limits[[id]]
        return(all(check_limit_value(setting, limit)))
        }
    )
    idx <- which(!check)
    if (length(idx) != 0) {
        msg <- paste("Whoops - A Limit was tighter than the constraint setting:",
                    paste(names(constraint.limits)[idx], collapse = ","))
        warning(msg)
        return(FALSE)
    } else {
        return(TRUE)
    }
}
#' Validity Check for DesignSettings.
#'
#' Validates whether a DesignSettings object has the correct structure.
#'
#' @param object A DesignSettings object to be checked for validity.
#' @return \code{TRUE} if \code{object} is valid, FALSE otherwise.
#' @keywords internal
check_settings_validity <- function(object) {
    # ensure that constraint limits are at least as wide as the constraint settings
    # check_limits shouldn't produce an error, just warn :-)
    check.limits <- check_limits(constraints(object@Input_Constraints), constraints(object@Input_Constraint_Boundaries))
	return(TRUE)
}
##########
# GETTERS
############

#' Getter/Setter for Constraints.
#'
#' Gets the active constraints of the provided
#' \code{DesignSettings} object \code{x}.
#'
#' @name constraints
#' @rdname constraints-methods
#' @exportMethod constraints
#' @keywords Settings
#' @family settings functions
setGeneric("constraints", function(x) standardGeneric("constraints"))

#' @rdname constraints-methods
#' @aliases constraints,DesignSettings-method
setMethod("constraints", "DesignSettings", function(x) {
    sel <- names(constraints(x@Input_Constraints))
    #print("Constraint getter:")
    #print(class(x@Input_Constraints))
    constraints(x@Input_Constraints)[con_select(sel)]
})

#' @rdname constraints-methods
#' @aliases constraints,DesignSettings-method
setMethod("constraints", c("AbstractConstraintSettings"), 
    function(x) {
        return(x@settings)    
    }
)

#' Getter/Setter for Coverage Constraints.
#'
#' Gets the coverage constraints of the provided
#' \code{DesignSettings} object \code{x}.
#' @name cvg_constraints
#' @exportMethod cvg_constraints
#' @rdname cvg_constraints-methods
#' @keywords Settings
#' @family settings functions
setGeneric("cvg_constraints", function(x) standardGeneric("cvg_constraints"))

#' @rdname cvg_constraints-methods
#' @aliases cvg_constraints,DesignSettings-method
setMethod("cvg_constraints", "DesignSettings", function(x) {
    sel <- names(constraints(x@Coverage_Constraints))
    # select only the possible constraints from the settings (tool dependencies):
    constraints(x@Coverage_Constraints)[select.constraints(sel)]
})

#' Getter for Filtering Constraints.
#'
#' Gets the constraints on the physicochemical properties
#' that are used for the filtering procedure when designing primers
#' using the \code{Input_Constraints} slot of the provided
#' \code{DesignSettings} object \code{x}.
#'
#' @name filters
#' @keywords internal
#' @rdname filters-methods
setGeneric("filters", function(x) standardGeneric("filters"))

#' @rdname filters-methods
#' @aliases filters,DesignSettings-method
#' @param x A \code{DesignSettings} object.
setMethod("filters", "DesignSettings", function(x) {
    # return all filters except for those only relevant to the optimization procedure:
    opti.names <- c("melting_temp_diff", "cross_dimerization")
    sel <- setdiff(names(constraints(x@Input_Constraints)), opti.names)
    constraints(x@Input_Constraints)[con_select(sel)]
})

#' Getter for Filtering Constraint Limits.
#'
#' Gets the limits on the constraints
#' that are used for the filtering procedure when designing primers
#' using the \code{Input_Constraint_Boundaries} slot of the provided
#' \code{DesignSettings} object \code{x}.
#'
#' @name filterLimits
#' @rdname filterLimits-methods
#' @keywords internal
setGeneric("filterLimits", function(x) standardGeneric("filterLimits"))

#' @rdname filterLimits-methods
#' @aliases filterLimits,DesignSettings-method
#' @param x A \code{DesignSettings} object.
setMethod("filterLimits", "DesignSettings", function(x) {
    opti.names <- c("melting_temp_diff", "cross_dimerization")
    sel <- setdiff(names(constraints(x@Input_Constraint_Boundaries)), opti.names)
    constraints(x@Input_Constraint_Boundaries)[con_select(sel)]
})

#' Getter for Optimization Constraints.
#'
#' Gets the constraints on the physicochemical properties
#' that are applied just before the optimization procedure
#' using the \code{Input_Constraints} slot of the provided
#' \code{DesignSettings} object \code{x}.
#'
#' @name opti
#' @rdname opti-methods
#' @keywords internal
setGeneric("opti", function(x) standardGeneric("opti"))

#' @rdname opti-methods
#' @aliases opti,DesignSettings-method
#' @param x A \code{DesignSettings} object.
setMethod("opti", "DesignSettings", function(x) {
    # define the constraints for filtering in the optimization procedure
    opti.names <- c("melting_temp_diff", "cross_dimerization")
    sel.constraints <- con_select(names(constraints(x@Input_Constraints)))[names(constraints(x@Input_Constraints)) %in% opti.names]
    opti.cons <- constraints(x@Input_Constraints)[sel.constraints]
    opti.cons
})

#' Getter for Optimization Constraint Limits.
#'
#' Gets the limits for the constraints
#' that are applied just before the optimization procedure
#' using the \code{Input_Constraint_Boundaries} slot of the provided
#' \code{DesignSettings} object \code{x}.
#'
#' @name optiLimits
#' @rdname optiLimits-methods
#' @keywords internal
setGeneric("optiLimits", function(x) standardGeneric("optiLimits"))

#' @rdname optiLimits-methods
#' @aliases optiLimits,DesignSettings-method
#' @param x A \code{DesignSettings} object.
setMethod("optiLimits", "DesignSettings", function(x) {
    opti.names <- c("melting_temp_diff", "cross_dimerization")
    sel <- names(constraints(x@Input_Constraint_Boundaries))[names(constraints(x@Input_Constraint_Boundaries)) %in% opti.names]
    opti.cons <- constraints(x@Input_Constraint_Boundaries)[sel]
    opti.cons
})

#' Getter/Setter for the PCR Conditions.
#'
#' Gets the PCR conditions that are defined in the provided
#' \code{DesignSettings} object \code{x}.
#'
#' @name PCR
#' @rdname PCR-methods
#' @exportMethod PCR
#' @family settings functions
#' @keywords Settings
setGeneric("PCR", function(x) standardGeneric("PCR"))

#' @rdname PCR-methods
#' @aliases PCR,DesignSettings-method
setMethod("PCR", "DesignSettings", function(x) {
    constraints(x@PCR_conditions)
})

#' Getter/Setter for Constraint Options.
#'
#' Gets the constraint settings of the provided 
#' \code{DesignSettings} object \code{x}.
#'
#' @name conOptions
#' @rdname conOptions-methods
#' @exportMethod conOptions
#' @keywords Settings
#' @family settings functions
setGeneric("conOptions", function(x) standardGeneric("conOptions"))

#' @rdname conOptions-methods
#' @aliases conOptions,DesignSettings-method
setMethod("conOptions", "DesignSettings", 
    function(x) {
        constraints(x@constraint_settings)
    }
)

#' Getter/Setter for Constraint Limits.
#'
#' Gets the constraint limits that are defined in the provided 
#' \code{DesignSettings} object \code{x}.
#' 
#' @name constraintLimits
#' @rdname constraintLimits-methods
#' @keywords Settings
#' @exportMethod constraintLimits
#' @family settings functions
setGeneric("constraintLimits", function(x) standardGeneric("constraintLimits"))

#' @rdname constraintLimits-methods
#' @aliases constraintLimits,DesignSettings-method
setMethod("constraintLimits", "DesignSettings", 
	function(x) {
        sel <- con_select(names(constraints(x@Input_Constraint_Boundaries)))
		constraints(x@Input_Constraint_Boundaries)[sel]
	}
)

###############
# SETTERS
###############

#' Setter for Constraints.
#'
#' Sets the active constraints of the provided
#' \code{DesignSettings} object \code{x}.
#' 
#' For an overview of permissible constraints,
#' please consider the \code{\link{ConstraintSettings}} documentation.
#'
#' @name constraints<-
#' @rdname constraints-methods
#' @exportMethod constraints<-
#' @param x A \code{DesignSettings} object.
#' @param value A list with constraint settings. Each list entry
#' should have a permissible name and consist of at most two
#' values providing the minimal and/or maximal allowed values, which
#' have to be denominated via \code{min} and \code{max}. 
#' @examples
#' # Load some settings
#' data(Ippolito)
#' # View the active constraints
#' constraints(settings)
#' # Require a minimal GC clamp extent of 0
#' constraints(settings)$gc_clamp["min"] <- 0
#' # View available constraints
#' settings
setGeneric("constraints<-", function(x, value) standardGeneric("constraints<-"))

#' @rdname constraints-methods
#' @aliases constraints,DesignSettings-method
setReplaceMethod("constraints", "DesignSettings", 
	function(x, value) {
        # modify constraint limits if necessary
        # a) ensure that constraints are in the right order
        con.order <- getOption("openPrimeR.constraint_order")
        sel.constraints <- con.order[con.order %in% names(value)]
        if (length(sel.constraints) != length(value)) {
           warning("The following constraint identifiers are not valid and were ignored: ",
                    paste0(names(value)[!names(value) %in% con.order], sep = ","), ".\nThe following constraints are valid: ", paste0(con.order, collapse = ","))
        }
        value <- value[con.order[con.order %in% names(value)]]
        # b) ensure that limit and constraint entries are concordant
        fixed.limits <- fix_constraint_boundaries(value, constraintLimits(x))
        constraints(x@Input_Constraint_Boundaries) <- fixed.limits
		constraints(x@Input_Constraints) <- value
		validObject(x)
		x
	}
)
#' @rdname constraints-methods
#' @aliases constraints,DesignSettings-method
setReplaceMethod("constraints", c("AbstractConstraintSettings", "list"),
    function(x, value) {
        m <- match(names(value), names(x@status))
        if (any(is.na(m))) {
            warning("The following constraints were ignored due to invalidity: ",
                    paste0(names(value)[is.na(m)], sep = ","),
                    ".\nOnly the follwing options are valid: ",
                    paste0(names(x@status), collapse = ","))
        }
        value <- value[!is.na(m)]
        x@settings <- value
        # activate the input constraints, deactivate all others
        x@status[names(value)] <- TRUE
        deactivate.idx <- setdiff(seq_along(x@status), m[!is.na(m)])
        x@status[deactivate.idx] <- FALSE
        # check for correct list structure with the validity function:
        validObject(x)
        return(x)
    }
)


#' Setter for Coverage Constraints.
#'
#' Sets the coverage constraints of the provided
#' \code{DesignSettings} object \code{x}.
#' 
#' @name cvg_constraints<-
#' @rdname cvg_constraints-methods
#' @exportMethod cvg_constraints<-
#' @param x A \code{DesignSettings} object.
#' @param value A list with coverage constraints. Each
#' list entry must have a permissible name and
#' contain a numeric vector with at most two components
#' describing the minimal and/or maximal required values that
#' are to be indicated via \code{min} and \code{max}.
#' The permissible contraint identifiers are documented in the
#' \code{\link{CoverageConstraints}} class.
#' @examples
#' # Load some settings
#' data(Ippolito)
#' # View all active coverage constraints
#' cvg_constraints(settings)
#' # Increase the maximal false positive rate to increase the sensitiviity of coverage predictions
#' cvg_constraints(settings)$coverage_model <- c("max" = 0.1) 
#' # View available coverage constraints:
#' settings
setGeneric("cvg_constraints<-", function(x, value) standardGeneric("cvg_constraints<-"))

#' @rdname cvg_constraints-methods
#' @aliases cvg_constraints,DesignSettings-method
setReplaceMethod("cvg_constraints", "DesignSettings", 
	function(x, value) {
		constraints(x@Coverage_Constraints) <- value
		validObject(x)
		x
	}
)
#' Setter for Constraint Limits.
#'
#' Sets the constraint limits of the provided 
#' \code{DesignSettings} object \code{x}.
#'
#' @rdname constraintLimits-methods
#' @exportMethod constraintLimits<-
#' @param x A \code{DesignSettings} object whose constraint limits
#' are to be modified.
#' @param value A list with constraint boundaries. The permissible fields 
#' of the list are provided in \code{\link{ConstraintSettings}}.
#' @examples
#' # Load some settings
#' data(Ippolito)
#' # View the active constraint limits
#' constraintLimits(settings)
#' # Extend the GC relaxation limit
#' constraintLimits(settings)$gc_clamp <- c("min" = 0, "max" = 6)
#' # View available constraints
#' settings
setGeneric("constraintLimits<-", function(x, value) standardGeneric("constraintLimits<-"))

#' @rdname constraintLimits-methods
#' @aliases constraintLimits,DesignSettings-method
setReplaceMethod("constraintLimits", "DesignSettings", 
	function(x, value) {
        # modify the settings if necessary
        # a) ensure that constraints are in the right order
        con.order <- getOption("openPrimeR.constraint_order")
        value <- value[con.order[con.order %in% names(value)]]
        # b) ensure that constraint limits and settings are still compatible
        fixed.settings <- fix_constraint_boundaries(constraints(x), value, fix.limit = FALSE)
		constraints(x@Input_Constraint_Boundaries) <- value
        constraints(x@Input_Constraints) <- fixed.settings
		validObject(x)
		x
	}
)

#' Setter for PCR Conditions.
#'
#' Sets the PCR conditions that are defined in the provided
#' \code{DesignSettings} object \code{x}.
#'
#' @name PCR<-
#' @rdname PCR-methods
#' @exportMethod PCR<-
#' @param x A \code{DesignSettings} object.
#' @param value A named list providing PCR conditions 
#' The permissible fields of the list and their types
#' are documented in the \code{\link{PCR_Conditions}} class.
#' @examples
#' # Load some settings
#' data(Ippolito)
#' # View the active PCR conditions
#' PCR(settings)
#' # Evaluate primers with a fixed annealing temperature
#' PCR(settings)$annealing_temperature <- 50 # celsius
#' # View available PCR conditions
#' settings

setGeneric("PCR<-", function(x, value) standardGeneric("PCR<-"))

#' @rdname PCR-methods
#' @aliases PCR,DesignSettings-method
setReplaceMethod("PCR", "DesignSettings", 
	function(x, value) {
		constraints(x@PCR_conditions) <- value
		validObject(x)
		x
	}
)

#' Setter for Constraint Options.
#'
#' Sets the constraint settings of the provided
#' \code{DesignSettings} object \code{x}.
#'
#' @name conOptions<-
#' @rdname conOptions-methods
#' @exportMethod conOptions<-
#' @param x A \code{DesignSettings} object.
#' @param value A list with constraint options. The permissible
#' fields of the list and their types are documented in the 
#' \code{\link{ConstraintOptions}} class.
#' @examples
#' # Load some settings
#' data(Ippolito)
#' # View the active constraint options
#' conOptions(settings)
#' # Prevent mismatch binding events
#' conOptions(settings)$allowed_mismatches <- 0
#' # View available constraint options
#' settings
setGeneric("conOptions<-", function(x, value) standardGeneric("conOptions<-"))

#' @rdname conOptions-methods
#' @aliases conOptions,DesignSettings-method
setReplaceMethod("conOptions", "DesignSettings", 
	function(x, value) {
		constraints(x@constraint_settings) <- value
		validObject(x)
		x
	}
)
#################################

#' Parse XML Constraint Data.
#'
#' Parses the constraint settings contained in an XML object.
#'
#' @param xml_data XML object from a parsed XML file.
#' @return List with constraint settings.
#' @keywords internal
parse.constraints <- function(xml_data) {
    result <- vector("list", length(xml_data))
    for (i in seq_along(xml_data)) {
        # for every set of constraints (filtering/optimization etc.)
        data <- xml_data[[i]]
        names(result)[i] <- data$name
        con.idx <- which(names(data) == "constraint")  # the entries we are looking for 
        values <- vector("list", length(con.idx))  # one entry for every constraint (melting temp, gc content etc.)
        for (j in seq_along(con.idx)) {
            idx <- con.idx[j]
            cur.con <- data[[idx]]
            names(values)[j] <- cur.con$name
            vals <- data[[idx]]$values
            # check whether constraint is min/max or value
            is.min.max.con <- any(grepl("=", vals))
            if (is.min.max.con) {
                v <- sapply(unlist(strsplit(vals, split = ",")), function(x) strsplit(x, 
                  split = "="))
                con.type <- sapply(v, function(x) x[1])
                settings <- sapply(v, function(x) as.numeric(x[2]))  # min/max something else setting?
                names(settings) <- con.type
            } else {
                settings <- suppressWarnings(as.numeric(vals))
                if (is.na(settings)) {
                  # vals was string (active/inactive setting)
                  settings <- vals  # don't transform to numeric
                }
            }
            values[[j]] <- settings
        }
        result[[i]] <- values
    }
    return(result)
}

#' Loading of Analysis Settings.
#'
#' Loads primer analysis settings from an XML file.
#'
#' If \code{filename} is not provided,
#' a default XMl settings file is loaded. Please review the 
#' function's examples to learn more about the default settings. If you want
#' to load custom settings, you can store a modified \code{DesignSettings}
#' object as an XML file using \code{\link{write_settings}}.
#'
#' @param filename Path to a valid XML file containing the 
#' primer analysis settings. By default, \code{filename} is set
#' to all settings that are shipped with openPrimeR and the lexicographically
#' first file is loaded.
#' @param frontend Indicates whether settings shall be loaded for the Shiny frontend. 
#' In this case no unit conversions for the PCR settings are performed.
#' The default setting is \code{FALSE} such that the correct units are used.
#' @return An object of class \code{DesignSettings}.
#' @family settings functions
#' @export
#' @examples
#' # Select the available settings
#' #' the available supplied settings by calling  
#' available.settings <- list.files(
#'      system.file("extdata", "settings", package = "openPrimeR"), 
#'      pattern = "*.xml", full.names = TRUE)
#' # Select one of the settings and load them
#' filename <- available.settings[1]
#' settings <- read_settings(filename)
#' # Modify, store, and read a settings object:
#' constraints(settings)$gc_clamp <- c("min" = 0, "max" = 5)
#' out.file <- tempfile("my_settings", fileext = ".xml")
#' write_settings(settings, out.file)
#' my_settings <- read_settings(out.file)
read_settings <- function(filename = list.files(
                             system.file("extdata", "settings", package = "openPrimeR"), 
                             pattern = "*.xml", full.names = TRUE),
                          frontend = FALSE) {
    # read xml constraint data from filename to constraint.settings objects
    if (length(filename) == 0) {
        stop("Please provide a non-empty filename.")
    }
    filename <- filename[1]
    if (!file.exists(filename)) {
        stop("Could not find settings file: '", filename, "'")
    }
    message("Reading settings file: ", basename(filename))
    # Validate XML prior to reading the settings
    schema.file <- system.file("extdata", "settings", 
                            "settings_schema.xsd", package = "openPrimeR")
    con.setting <- try({
        xsd <- XML::xmlTreeParse(schema.file,
                                isSchema = TRUE, useInternal = TRUE)
        doc <- XML::xmlInternalTreeParse(filename)
        xml.check <- XML::xmlSchemaValidate(xsd, doc)
        if (xml.check$status != 0) {
            msg <- xml.check$errors
            stop(paste("The settings XML file ", filename, 
                 " did not adhere to the required schema. The error was: ",
                 msg, sep = ""))
        }
        data <- XML::xmlParse(filename)
        xml_data <- XML::xmlToList(data)
        con.setting <- parse.constraints(xml_data)
        con.setting
    })
    if (class(con.setting) == "try-error") {
        msg <- "Error while parsing settings XML file. Please check your input."
        my.error("XML_Parsing_Error", msg)
    }
    # convert active/inactive to boolean
    # for each slot, convert all 'active'/'inactive' annotations to TRUE/FALSE
    for (i in seq_along(con.setting)) {
        c.setting <- lapply(con.setting[[i]], function(x)  {
            if (length(x) == 1 && x[[1]] %in% c("active", "inactive")) {
                return(ifelse(x[[1]] == "active", TRUE , ifelse(x[[1]] == "inactive", FALSE, x)))
            } else {
                return(x)
            }
        })
        con.setting[[i]] <- c.setting
    }
    if (!frontend) {
        # for backend: convert the settings to non-prefixed notation
        con.setting$PCR_conditions <- convert.PCR.units(con.setting$PCR_conditions, to.mol = TRUE)
    }
	# construct DesignSettings object
    settings <- DesignSettings(con.setting$Input_Constraints,
               con.setting$Input_Constraint_Boundaries, 
               con.setting$Coverage_Constraints, 
               con.setting$PCR_conditions, con.setting$constraint_settings)
    return(settings)
}
#' Conversion of PCR Units
#'
#' Converts frontend PCR concentration units to the units used for the backend.
#'
#' @param pcr.settings List with several PCR settings (concentrations).
#' @param to.mol If \code{TRUE}, convert to the molar concentration.
#' If \code{FALSE} convert to the unit representation in the XML.
#' @return List with concentrations for usage in the backend.
#' @keywords internal
convert.PCR.units <- function(pcr.settings, to.mol = TRUE) {
    # Ion concentrations are in mM in XML
    ions <- c("Na_concentration", "Mg_concentration", "K_concentration",
            "Tris_concentration")
    if (!all(ions %in% names(pcr.settings))) {
        stop("Concentration names have changed in xml.")
    }
    for (i in seq_along(ions)){
        if (to.mol) {
            # convert from mM to molar
            pcr.settings[[ions[i]]] <- pcr.settings[[ions[i]]] * 1e-3
        } else {
            # convert from molar to mM
            pcr.settings[[ions[i]]] <- pcr.settings[[ions[i]]] / 1e-3
        }
    }
    # Primer/template concentrations are in nM in XML
    other <- c("primer_concentration", "template_concentration")
    if (!all(other %in% names(pcr.settings))) {
        stop("Concentration names have changed in xml.")
    }
    for (i in seq_along(other)){
        if (to.mol) {
            pcr.settings[[other[i]]] <- pcr.settings[[other[i]]] * 1e-9
        } else {
            pcr.settings[[other[i]]] <- pcr.settings[[other[i]]] / 1e-9
        }
    }
    return(pcr.settings)
}

#' Gather all Coverage Constraints.
#'
#' Constructor for coverage constraint settings.
#'
#' @param allowed.stop.codons Whether mismatch binding events inducing
#' stop codons in the amino acid sequence are allowed.
#' @param allowed.efficiency Min/max for primer efficiency.
#' @param disallowed.mismatch.pos The positions from the 3' terminal
#' end of primers where mismatches shall be prevented.
#' @param allowed.anneal.deltaG Maximal allowed free energy of template-primer annealing.
#' @param allowed.substitutions Whether mismatch binding events inducing substitutions
#' in the amino acid sequence are allowed.
#' @return List with all coverage constraint settings.
#' @keywords internal
get.cvg.constraint.settings <- function(allowed.stop.codons, allowed.efficiency, 
                                        disallowed.mismatch.pos, allowed.anneal.deltaG,
                                        allowed.substitutions, allowed.coverage.model) {

    if (allowed.stop.codons) {
        # stop codons allowed
        #allowed.stop.range <- c("min" = 0, "max" = 1)
        allowed.stop.range <- NULL # don't consider constraint
    } else {
        # no stop codons allowed
        allowed.stop.range <- c("min" = 0, "max" = 0) 
    }
    if (allowed.substitutions) {
        # substitutions allowed
        #allowed.sub.range <- c("min" = 0, "max" = 1)
        allowed.sub.range <- NULL # don't consider constraint
    } else {
        # no substitutions allowed
        allowed.sub.range <- c("min" = 0, "max" = 0) 
    }
    settings <- list()
    if (length(allowed.stop.codons) != 0 && length(allowed.stop.range) != 0) {
        settings$stop_codon <- allowed.stop.range
    }
    if (length(allowed.efficiency) != 0) {
        settings$primer_efficiency <- allowed.efficiency
    }
    if (length(allowed.anneal.deltaG) != 0) {
        settings$annealing_DeltaG <- allowed.anneal.deltaG
    }
    if (length(disallowed.mismatch.pos) != 0) {
        settings$terminal_mismatch_pos <- c(disallowed.mismatch.pos + 1)
    }
    if (length(allowed.sub.range) != 0 && length(allowed.sub.range) != 0) {
        settings$substitution <- allowed.sub.range
    }
    if (length(allowed.coverage.model) != 0) {
        settings$coverage_model <- allowed.coverage.model
    }
    return(settings)
} 

#' Gather all Other Constraints (for Shiny frontend).
#'
#' Constructor for other constraint settings (non-PCR, non-filtering, non-optimization).
#'
#' @param allowed_mismatches Allowed mismatches for primers binding events.
#' @param allowed_other_binding_ratio Ratio of primers allowed to bind to non-target regions.
#' @param allowed_region_definition The definition of the allowed region.
#' @return List with all other constraint settings.
#' @keywords internal
get.other.constraint.settings <- function(allowed_mismatches, 
    allowed_other_binding_ratio, allowed_region_definition) {
  
    settings <- list(allowed_mismatches = allowed_mismatches, 
        allowed_other_binding_ratio = allowed_other_binding_ratio, 
        allowed_region_definition = allowed_region_definition)
    #if (hexamer_coverage == "active") { 
        ## don't add hexamer coverage as constraint if inactive anyway
        #settings <- settings["hexamer_coverage"] <- "active"
    #}
    return(settings)
} 
#' Constraint XML Format.
#'
#' Format constraint settings for XML output.
#'
#' @param constraints List with constraint settings.
#' @param set.name Identifier for the constraint settings.
#' @return XML string containing the constraint settings.
#' @keywords internal
constraints.xml.format <- function(constraints, set.name) {
    # create a consistent format for the constraints for xml output
    main <- "constraint"
    out <- rep("", length(constraints))
    for (i in seq_along(constraints)) {
        # message(constraints[[i]])
        if (length(names(constraints[[i]])) != 0) {
            val.string <- paste(names(constraints[[i]]), "=", constraints[[i]], collapse = ",", 
                sep = "")
        } else {
            # no min/max -> just the value itself
            val.string <- paste(constraints[[i]], collapse = ",", sep = "")
        }
        # warning: saveXML/listToXML creates many textConnections and doesn't close them!!
        # realiziation: I don't even need the XML::saveXML call after using listToXml
        a <- xmlToChar(listToXml(names(constraints)[[i]], "name"))
        b <- xmlToChar(listToXml(val.string, "values"))
        out[i] <- paste(a, b, sep = "\n")
    }
    # force closing of all connections: this could be very bad if any other connection is open for a reason...
    if (length(out) != 0) {
        result <- paste("<constraint>", out, "</constraint>", collapse = "\n", sep = "\n")
    } else {
        result <- ""
    }
    result <- paste("<constraint_set>\n", "<name>", set.name, "</name>\n", result, 
        "</constraint_set>", sep = "")
    return(result)
}

#' Storing Design Settings to Disk.
#'
#' Stores primer analysis settings to a file in XML format.
#'
#' @param settings A \code{DesignSettings} object to be stored to disk.
#' @param filename A character vector specifying the location 
#' where the settings should be stored as an XML file.
#' @return Outputs the return status from closing the connection.
#' @export
#' @keywords Settings
#' @examples
#' xml <- settings.xml <- system.file("extdata", "settings", 
#'        "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
#' settings <- read_settings(xml)
#' out.file <- tempfile("my_settings", fileext = ".xml")
#' write_settings(settings, out.file)
write_settings <- function(settings, filename) {
    XML <- create.constraint.XML(constraints(settings), constraintLimits(settings), 
                                 cvg_constraints(settings), PCR(settings), conOptions(settings))
    f <- file(filename)
    write(XML, f)
    return(close(f))
}
#' XML Output of Constraints
#'
#' Creates an XML summarizing all settings.
#'
#' @param filtering.constraints List with constraint settings for filtering.
#' @param c.f.lim Relaxation limits for the filtering constraints.
#' @param cvg.constraints List with constraints for coverage computations.
#' @param PCR.settings Settings for the PCR.
#' @param constraint.settings Other settings of constraints (e.g. coverage).
#' @return String in XML format containing all constraint settings.
#' @keywords internal
create.constraint.XML <- function(filtering.constraints, c.f.lim, 
                                  cvg.constraints, PCR.settings, constraint.settings) {
    comment <- "<!--These settings were automatically generated by openPrimeR. Please do not edit!-->"
    root.node <- "<design_settings>"
    end.node <- "</design_settings>"
    prefix <- character(0)
    ###############
    # Constraints
    ##############
    if (length(filtering.constraints) != 0) {
        m <- match(names(filtering.constraints), names(c.f.lim))
        c.f.lim <- c.f.lim[m[!is.na(m)]]
    }
    if (length(filtering.constraints) != 0) {
        c.f.in <- constraints.xml.format(filtering.constraints, "Input_Constraints")
    } else {
        c.f.in <- ""
    }
    ######################
    # Constraints limits
    ######################
    if (length(c.f.lim) != 0) {
        c.f.lim <- constraints.xml.format(c.f.lim, "Input_Constraint_Boundaries")
    } else {
        c.f.lim <- ""
    }
    ######################
    # Coverage constraints
    ######################
    if (length(cvg.constraints) != 0) {
        c.cvg <- constraints.xml.format(cvg.constraints, "Coverage_Constraints")
    } else {
        c.cvg <- ""
    }
    #############
    # PCR settings
    ##############
    # convert the raw, molar concentrations to XML units (nM and such)
    PCR.settings <- convert.PCR.units(PCR.settings, to.mol = FALSE)
    # turn TRUE/FALSE to active/inactive
    PCR.settings <- lapply(PCR.settings, function(x) ifelse(is.logical(x[[1]]), 
                                  ifelse(x[[1]] == TRUE, "active", ifelse(x[[1]] == FALSE, "inactive", x[[1]])), 
                                  x[[1]]))
    PCR.settings <- constraints.xml.format(PCR.settings, "PCR_conditions")
    #############################
    # Constraint options
    #######################
    # need to replace TRUE/FALSE with active/inactive
    constraint.settings <- lapply(constraint.settings, function(x) ifelse(is.logical(x[[1]]), 
                                  ifelse(x[[1]] == TRUE, "active", ifelse(x[[1]] == FALSE, "inactive", x[[1]])), 
                                  x[[1]]))
    constraint.settings <- constraints.xml.format(constraint.settings, "constraint_settings")  # constraint options
    sep <- "\n"
    out <- paste(comment, root.node, c.f.in, c.f.lim, c.cvg, PCR.settings, 
                 constraint.settings, end.node, sep = sep)
    return(out)
}

#' Conversion of XML to Character.
#'
#' Converts an XML object to a character string.
#'
#' @param xml An xml object to be converted to character.
#' @return A character vector.
#' @keywords internal
xmlToChar <- function(xml) {
    file <- textConnection(NULL, "w")
    sink(file)
    on.exit({sink();close(file)})
    # we need this print statement here:
    print(xml)
    paste(textConnectionValue(file), collapse = "\n")
}
#' Convert List to XML.
#'
#' Can convert list or other object to an xml object using xmlNode.
#'
#' @title List to XML
#' @param item 
#' @param tag xml tag
#' @return xmlNode
#' @author David LeBauer, Carl Davidson, Rob Kooper
#' @keywords internal
listToXml <- function(item, tag) {
    # just textnode, or empty node with attributes
    if (typeof(item) != "list") {
        if (length(item) > 1) {
            xml <- XML::xmlNode(tag)
            for (name in names(item)) {
                XML::xmlAttrs(xml)[[name]] <- item[[name]]
            }
            return(xml)
        } else {
            return(XML::xmlNode(tag, item))
        }
    }
    
    # create the node
    if (identical(names(item), c("text", ".attrs"))) {
        # special case a node with text and attributes
        xml <- XML::xmlNode(tag, item[["text"]])
    } else {
        # node with child nodes
        xml <- XML::xmlNode(tag)
        for (i in 1:length(item)) {
            if (names(item)[i] != ".attrs") {
                xml <- XML::append.xmlNode(xml, listToXml(item[[i]], names(item)[i]))
            }
        }
    }
    
    # add attributes to node
    attrs <- item[[".attrs"]]
    for (name in names(attrs)) {
        XML::xmlAttrs(xml)[[name]] <- attrs[[name]]
    }
    return(xml)
}
#' Gather all PCR settings.
#'
#' Gathers all PCR settings (e.g. for XML output).
#'
#' @param annealing_temp Annealing temperature in Celsius.
#' @param Na_concentration Sodium ion concentration.
#' @param Mg_concentration Magensium ion concentration.
#' @param K_concentration Potassium ion concentration.
#' @param Tris_concentration Tris ion concentration.
#' @param primer_concentration Primer concentration.
#' @param template_concentration Template concentration.
#' @return List with all PCR settings.
#' @keywords internal
get.PCR.settings <- function(use_taq_polymerase, annealing_temp, Na_concentration, Mg_concentration, 
    K_concentration, Tris_concentration, primer_concentration, template_concentration, nbr.cycles) {
        settings <- list(use_taq_polymerase = use_taq_polymerase, Na_concentration = Na_concentration, 
        Mg_concentration = Mg_concentration, K_concentration = K_concentration, Tris_concentration = Tris_concentration, 
        primer_concentration = primer_concentration, template_concentration = template_concentration, cycles = nbr.cycles)
    if (length(annealing_temp) != 0) {
        settings <- c("annealing_temp" = annealing_temp, settings)
    }
    return(settings)
}

#' Conversion of Constraints List to Data Frame.
#'
#' Converts the input constraints to a data frame representation.
#'
#' @param limit.constraints A list with constraints.
#' @param out.names The desired column names.
#' @param format.type The type of formatting to be performed on the table
#' @return A data frame giving an overview of the constraints.
#' @keywords internal
constraints.to.df <- function(limit.constraints, out.names, format.type = c("backend", "shiny", "report")) {
    if (length(out.names) != 1) {
        stop("out.names must have length 1")
    }
    format.type <- match.arg(format.type)
    cnames <- c("Constraint", "min", "max")
    df <- NULL
    for (i in seq_along(limit.constraints)) {
        c <- limit.constraints[[i]]
        con.name <- names(limit.constraints)[i]
        criteria <- names(c)
        idx <- sapply(criteria, function(x) which(x == cnames))
        idx.length <- unlist(lapply(idx, length))
        sel <- which(idx.length != 0)  # add sel to df
        n.sel <- which(idx.length == 0)
        nrows <- max(idx.length) + length(n.sel)
        d <- data.frame(matrix(rep(NA, nrows * length(cnames)), nrow = nrows))
        colnames(d) <- cnames
        rnames <- rep(NA, nrow(d))
        if (length(sel) != 0) {
            rnames[1] <- con.name
            dir <- criteria[sel]
            d[1, dir] <- c[sel]
        }
        if (length(n.sel) != 0) {
            origin.name <- criteria[n.sel]
            dir <- ifelse(grepl("min", origin.name), "min", "max")
            name <- sapply(seq_along(origin.name), function(x) gsub(paste(".", dir[x], 
                sep = ""), "", paste(con.name, ":", origin.name[x], sep = "")))
            rnames[n.sel] <- name
            d.s <- ifelse(length(sel) != 0, 2, 1)
            d.e <- d.s + length(n.sel) - 1
            for (j in seq_along(dir)) {
                d[(d.s:d.e)[j], dir[j]] <- c[n.sel][j]
            }
        }
        rownames(d) <- rnames
        df <- rbind(df, d)
    }
    if (length(df) != 0) {
        df$Constraint <- rownames(df)
        if (format.type == "shiny") {
            # annotate units with html units
            df$Constraint <- unlist(constraints_to_unit(df$Constraint, TRUE, "HTML"))
        } else if (format.type == "report") {
            df$Constraint <- unlist(constraints_to_unit(df$Constraint, TRUE, "report"))
        } else {
            # don't annotate units
            #df$Constraint <- unlist(constraints_to_unit(df$Constraint, FALSE)) # keep the input fields
        }
        df[,2:3] <- apply(df[,2:3], 2, function(x) round(x,3))
        # introduce interval notation:
        if (format.type %in% c("shiny", "report")) {
            # change NA to infinity symbols
            if (format.type == "shiny") {
                s1 <- "-&infin;"
                s2 <- "&infin;"
            } else {
                s1 <- "$-\\infty$"
                s2 <- "$\\infty$"
            }
            df[,2][is.na(df[,2])] <- s1
            df[,3][is.na(df[,3])] <- s2
        }
        df[,2] <- paste0("[", df[,2], ", ", df[,3], "]")
        df <- df[, -3] # remove the third column
        # modify column name of 2nd column
        colnames(df)[2] <- out.names
    }
    rownames(df) <- NULL
    return(df)
}

#' Renaming of Constraint Options.
#'
#' Renames the input list with constraint options.
#'
#' @param constraint.options A list with constraint options.
#' @return A list with renamed constraint options.
#' @keywords internal
rename.constraint.options <- function(constraint.options) {
    new.names <- Hmisc::capitalize(gsub("_", " ", names(constraint.options)))
    names(constraint.options) <- new.names
    return(constraint.options)
}
#' Output a Constraint Overview Table
#'
#' Outputs a table showing the values of constraints.
#'
#' @param constraints List with constraint settings.
#' @param constraint.limits List with constraint limits.
#' @param constraints.used.fw Constraints used for forward primer design.
#' @param constraints.used.rev Constraints used for reverse primer design.
#' @param format.type The type of formatting to be performed on entries.
#' @return Data frame with summary of constraints.
#' @keywords internal
create.constraint.table <- function(constraints, constraint.limits = NULL, constraints.used.fw = NULL, 
                                    constraints.used.rev = NULL, format.type = c("backend", "shiny", "report")) {
    if (length(constraints) == 0) {
        return(NULL)
    }
    format.type <- match.arg(format.type)
    input.constraints <- constraints
    limit.constraints <- constraint.limits
    output.constraints.fw <- constraints.used.fw
    output.constraints.rev <- constraints.used.rev
    if (length(input.constraints) != 0) {
        input.constraints <- constraints.to.df(input.constraints, "Target range", format.type)
    }
    if (length(output.constraints.fw) != 0) {
        out.names <- "Used range (fw)"
        # reorder
        m <- match(names(constraints), names(output.constraints.fw))
        output.constraints.fw <- constraints.to.df(output.constraints.fw[m[!is.na(m)]], out.names, format.type)
    }
    if (length(output.constraints.rev) != 0) {
        out.names <- "Used range (rev)"
        # reorder
        m <- match(names(constraints), names(output.constraints.rev))
        output.constraints.rev <- constraints.to.df(output.constraints.rev[m[!is.na(m)]], out.names, format.type)
    }
    if (length(input.constraints) != 0 && length(limit.constraints) != 0) {
        # only show relevant limits/reorder
        m <- match(names(constraints), names(constraint.limits))
        limit.constraints <- limit.constraints[m[!is.na(m)]]
        out.names <- "Limit range"
        limit.constraints <- constraints.to.df(limit.constraints, out.names, format.type)
    }
    result <- input.constraints
    if (length(limit.constraints) != 0) {
        result <- cbind(result, limit.constraints[,which(colnames(limit.constraints) != "Constraint"), drop = FALSE])
    }
    if (length(output.constraints.fw) != 0) {
        result <- cbind(result, output.constraints.fw[, which(colnames(output.constraints.fw) != "Constraint"), drop = FALSE])
    }
    if (length(output.constraints.rev) != 0) {
        result <- cbind(result, output.constraints.rev[, which(colnames(output.constraints.rev) != "Constraint"), drop = FALSE])
    }
    if (length(result) != 0) {
        # reorder: put Constraint name first
        o <- which(colnames(result) == "Constraint")
        o <- c(o, which(colnames(result) != "Constraint"))
        result <- result[,o]
    }
    rownames(result) <- NULL
    return(result)
}

#' Creation of a Table for Other Constraint Settings.
#'
#' @param other.settings List with other constraint settings.
#' @param format.type How the table shall be formatted.
#' @return A data frame.
#' @keywords internal
create.other.table <- function(other.settings, col.names, format.type) {
    option.names <- names(other.settings)
    if (format.type == "backend") {
        # nothing to format
    } else if (format.type == "shiny") {
        option.names <- constraints_to_unit(option.names, TRUE, "HTML")
    } else { # report format
        option.names <- constraints_to_unit(option.names, TRUE, format.type)
    }
    option.names <- unname(unlist(option.names))
    df <- data.frame(option.names, unlist(as.character(other.settings)))
    if (length(col.names) != ncol(df)) {
        stop("Dimensions of col.names and df do not agree.")
    }
    colnames(df) <- col.names
    return(df)
}
#' Creation of a Table for Constraint Options.
#'
#' @param other.settings List with constraint options
#' @param format.type How the table shall be formatted.
#' @return A data frame.
#' @keywords internal
create.options.table <- function(other.settings, 
                        format.type = c("backend", "shiny", "report")) {
    format.type <- match.arg(format.type)
    df <- create.other.table(other.settings, c("Option", "Setting"), format.type)
    rownames(df) <- NULL
    return(df)
}
#' Creation of a Table for PCR Conditions
#'
#' @param other.settings List with PCR settings.
#' @param format.tyep How the table shall be formatted.
#' @return A data frame.
#' @keywords internal
create.PCR.table <- function(other.settings, 
                    format.type = c("backend", "shiny", "report")) {
    format.type <- match.arg(format.type)
    df <- create.other.table(other.settings, c("Condition", "Setting"), format.type)

    rownames(df) <- NULL
    return(df)
}
#' Constraint list comparison
#'
#' Determines whether two list with constraints are identical.
#'
#' @param A First constraint list.
#' @param B Second constraint list.
#'
#' @return TRUE if the constraints are identical, FALSE else.
#' @keywords internal
compare.constraints <- function(A, B) {
    if (length(A) != length(B)) {
        return(FALSE)
    }
    m <- match(names(A), names(B))
    if (any(is.na(m))) {
        # names don't match
        return(FALSE)
    }
    # need to check further because the names of all lists agree
    for (i in seq_along(A)) {
        if (any(A[[i]] != B[[m[i]]])) {
            return(FALSE)
        }
    }
    return(TRUE)
}
#' Retrieval of Tool Information.
#'
#' Constructs a data frame containing information about the tools.
#'
#' @return A data frame with information about the required tools.
#' @keywords internal
get.static.tool.info <- function() {
    tools <- c("MELTING", "ViennaRNA", "OligoArrayAux", "MAFFT",
              "Selenium", "Pandoc", "PhantomJS")
    blank <- rep(NA, length(tools))
    names(blank) <- tools
    purposes <- blank
    purposes["MELTING"] <- "Melting temperatures"
    purposes["ViennaRNA"] <- "Secondary structures"
    purposes["OligoArrayAux"] <- "Primer efficiencies and dimerization"
    purposes["MAFFT"] <- "Multiple sequence alignments"
    purposes["Selenium"] <- "IMGT queries"
    purposes["Pandoc"] <- "PDF reports"
    purposes["PhantomJS"] <- "IMGT queries"
    # complete purposes, TODO
    tool.ex <- blank
    tool.ex["MELTING"] <- "melting-batch"
    tool.ex["ViennaRNA"] <- "RNAfold"
    tool.ex["OligoArrayAux"] <- "hybrid-min"
    tool.ex["MAFFT"] <- "mafft"
    tool.ex["Selenium"] <- "Python"
    tool.ex["Pandoc"] <- "Pandoc"
    tool.ex["PhantomJS"] <- "phantomjs"
    locations <- sapply(tool.ex, Sys.which)
    tool.URLs <- blank
    tool.URLs["MELTING"] <- "http://www.ebi.ac.uk/biomodels/tools/melting/"
    tool.URLs["ViennaRNA"] <- "http://www.tbi.univie.ac.at/RNA/"
    tool.URLs["OligoArrayAux"] <- "http://unafold.rna.albany.edu/OligoArrayAux.php"
    tool.URLs["MAFFT"] <- "http://mafft.cbrc.jp/alignment/software/"
    tool.URLs["Selenium"] <- "http://selenium-python.readthedocs.io/"
    tool.URLs["Pandoc"] <- "http://pandoc.org"
    tool.URLs["PhantomJS"] <- "http://phantomjs.org/"
    result <- data.frame("Tool" = names(purposes), "Purpose" = purposes, "Executable" = tool.ex, "URL" = tool.URLs)
    return(result)
}
#' Creation of an Overview of Third-Party Tools.
#'
#' Creates a table of required third-party tools and their installation status.
#'
#' @param AVAILBLE.TOOLS A vector whose names give the required tools and whose
#' entries give their installation status as logicals.
#' @param If \code{for.shiny} is \code{TRUE}, provide the URLs for the tool using HTML.
#' @return A data frame with information on third-part tools.
#' @keywords internal
build.tool.overview <- function(AVAILABLE.TOOLS, for.shiny = FALSE) {
    status <- ifelse(AVAILABLE.TOOLS, "Available", "Unavailable")
    tool.info <- get.static.tool.info()
    m <- match(tool.info$Tool, names(status))
    status <- status[m]
    tool.df <- data.frame(Tool = tool.info$Tool, Status = status, Purpose = tool.info$Purpose, 
                          Executable = tool.info$Executable, URL = tool.info$URL, stringsAsFactors = FALSE)
    if (for.shiny) {
        # add URLs to names
        tool.names <- paste0("<a href='", tool.info$URL, "' target='_blank'>", tool.info$URL, "</a>")
        tool.df$Tool <- tool.names
        # remove URL column
        tool.df <- tool.df[,colnames(tool.df) != "URL"]
    }
    return(tool.df)
}
#' Selection of Constraints.
#'
#' Selects constraints that can be computed according to installed third-party software.
#' This function is only used for initializing the 'constraint_order' option.
#'
#' @param active.constraints A vector whose names give the constraints to be checked.
#' @return A vector of useable constraint identifiers.
#' @keywords internal
select.constraints <- function(active.constraints) {
    tool.info <- build.tool.overview(check.tool.function(), for.shiny = FALSE)
    # check for each tool:
    melting.available <- tool.info[tool.info$Tool == "MELTING", "Status"] == "Available"
    vienna.available <- tool.info[tool.info$Tool == "ViennaRNA", "Status"] == "Available"
    oligo.available <- tool.info[tool.info$Tool == "OligoArrayAux", "Status"] == "Available"
    rm.constraints <- NULL
    #melting.constraints <- c("melting_temp_range", "melting_temp_diff")
    melting.constraints <- NULL # melting temp is now computed by empiric formula if melting is not present
    # this ensures we can also compute an annealing temperature (required for some constraints ..)
    vienna.constraints <- c("secondary_structure")
    oligo.constraints <- c("primer_efficiency", "self_dimerization", "cross_dimerization", "annealing_DeltaG", "coverage_model")
    if (!melting.available) {
        rm.constraints <- c(rm.constraints, melting.constraints)
    }
    if (!vienna.available) {
        rm.constraints <- c(rm.constraints, vienna.constraints)
    }
    if (!oligo.available) {
        rm.constraints <- c(rm.constraints, oligo.constraints)
    }
    ignore.constraints <- intersect(active.constraints, rm.constraints)
    new.constraints <- setdiff(active.constraints, ignore.constraints)
    return(new.constraints)
}

#' Quick Selection of Constraints.
#'
#' Select constraints that can be used according to third-party tools quickly.
#'
#' @param active.constraints Identifiers of constraints.
#' @return The identifiers of constraints that can be computed.
#' @keywords internal
con_select <- function(active.constraints) {
    out <- active.constraints[active.constraints %in% getOption("openPrimeR.constraint_order")]
    return(out)
}

tool.info <- get.static.tool.info()
#' Format Constraint Names.
#'
#' Formats constraint names for frontend output.
#'
#' @param constraints The character vector of constraints to transform.
#' @return A character vector with formatted constraint names.
#' @keywords internal
format.constraints <- function(constraints) {
    mapping <- list("primer_coverage" = "Coverage", 
        "primer_specificity" = "Specificity",
        "primer_length" = "Length",
        "gc_clamp" = "GC clamp",
        "gc_ratio" = "GC ratio",
        "no_runs" = "Runs",
        "no_repeats" = "Repeats", 
        "self_dimerization" = "Self dimers",
        "cross_dimerization" = "Cross dimers",
        "melting_temp_range" = "Tm range",
        "secondary_structure" = "Structures",
        "melting_temp_diff" = "Tm deviation",
        # coverage constraints
        "primer_efficiency" = "Efficiency",
        "annealing_DeltaG" = "Annealing",
        "stop_codon" = "Stop codons",
        "substitution" = "Substitutions",
        "terminal_mismatch_pos" = "3' Mismatch Position",
        "coverage_model" = "Coverage Model FPR",
        # PCR settings
        "use_taq_polymerase" = "Taq polymerase",
        "Na_concentration" = "[Na]",
        "Mg_concentration" = "[Mg]",
        "K_concentration" = "[K]",
        "Tris_concentration" = "[Tris buffer]",
        "primer_concentration" = "[Primer]",
        "template_concentration" = "[Template]",
        "cycles" = "PCR cycles",
        "annealing_temp" = "Annealing temperature",
        # Other options
        "allowed_mismatches" = "Allowed mismatches",
        "allowed_other_binding_ratio" = "Allowed off-target binding ratio",
        "allowed_region_definition" = "Binding region definition")
    m <- match(constraints, names(mapping))
    idx <- which(is.na(m))
    out <- mapping[m]
    if (length(idx) != 0) {
        msg <- paste("Could not format the following constraints:",
                paste(constraints[idx], collapse = ","))
        warning(msg)
        # keep original name if not found
        out[idx] <- constraints[idx]
    }
    return(out)
}
#' Mapping of Constraints to Units.
#'
#' Maps constraints to units for plotting.
#'
#' @param constraint The names of the constraints to convert to their plot identifiers (units).
#' @param use.unit Whether constraint names should be annotated with their units.
#' @param use.HTML Whether constraint units should be annotated with HTML units.
#' @return A list of constraint names.
#' @keywords internal
constraints_to_unit <- function(constraint, use.unit = TRUE, 
                    format.type = c("backend", "HTML", "report")) {
    format.type <- match.arg(format.type)
    if (length(constraint) == 0) {
        # nothing to modify
        return(constraint)
    }
    mod.con <- format.constraints(constraint)
    unit.mapping.R <- list(
        melting_temp_range = expression(paste("T"[m], " [", degree * 
                                        C, "]", sep = "")),
        melting_temp_diff = expression(paste("T"[m], " deviation [", degree * 
                                        C, "]", sep = "")),
        self_dimerization = expression(paste("Self dimer ", Delta, " G [kcal/mol]", sep = "")),
        cross_dimerization = expression(paste("Cross dimer ", Delta, " G [kcal/mol]", sep = "")),
        secondary_structure = expression(paste("Structure ", Delta, " G [kcal/mol]", sep = "")),
        annealing_DeltaG = expression(paste("Annealing ", Delta, " G [kcal/mol]", sep = ""))
        )
    unit.mapping.report <- list(
        melting_temp_range = "T\\textsubscript{m} range [\\textdegree C]",
        melting_temp_diff = "T\\textsubscript{m}  deviation [\\textdegree C]",
        self_dimerization = "Self dimer $\\Delta$G[$\\frac{\\text{kcal}}{\\text{mol}}$]",
        cross_dimerization = "Cross dimer $\\Delta$G[$\\frac{\\text{kcal}}{\\text{mol}}$]",
        secondary_structure = "Structure $\\Delta$G[$\\frac{\\text{kcal}}{\\text{mol}}$]",
        annealing_DeltaG = "Annealing $\\Delta$G[$\\frac{\\text{kcal}}{\\text{mol}}$]",
        # PCR conditions
        Na_concentration = "[Na\\textsuperscript{+}] [M]",
        Mg_concentration = "[Mg\\textsuperscript{2+}] [M]",
        K_concentration = "[K\\textsuperscript{+}] [M]",
        Tris_concentration = "[Tris buffer] [M]",
        primer_concentration = "[Primer] [M]",
        template_concentration = "[Template] [M]",
        annealing_temp = "Annealing temperature [\\textdegree C]"
        )
    unit.mapping.html <- list(
        melting_temp_range = "T<sub>m</sub> range [&#x2103;]",
        melting_temp_diff = "T<sub>m</sub> deviation [&#x2103;]",
        self_dimerization = "Self dimer &Delta;G[<sup>kcal</sup>&frasl;<sub>mol</sub>]",
        cross_dimerization = "Cross dimer &Delta;G[<sup>kcal</sup>&frasl;<sub>mol</sub>]",
        secondary_structure = "Structure &Delta;G[<sup>kcal</sup>&frasl;<sub>mol</sub>]",
        annealing_DeltaG = "Annealing &Delta;G[<sup>kcal</sup>&frasl;<sub>mol</sub>]",
        # PCR conditions
        Na_concentration = "[Na<sup>+</sup>] [M]",
        Mg_concentration = "[Mg<sup>2+</sup> [M]",
        K_concentration = "[K<sup>+</sup>] [M]",
        Tris_concentration = "[Tris buffer] [M]",
        primer_concentration = "[Primer] [M]",
        template_concentration = "[Template] [M]",
        annealing_temp = "Annealing temperature [&#x2103;]"
        )
    if (use.unit) {
        if (format.type == "HTML") {
            mapping <- unit.mapping.html
        } else if (format.type == "report") {
            mapping <- unit.mapping.report
        } else {
            mapping <- unit.mapping.R
        }
        m <- match(constraint, names(mapping))
        idx <- which(!is.na(m))
        if (length(idx) != 0) {
            mod.con[idx] <- lapply(idx, function(x) mapping[[m[x]]])
        }
    }
    return(mod.con)
}


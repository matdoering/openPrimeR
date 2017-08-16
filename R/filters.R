#' Filter Multiple Primer Sets.
#' 
#' Filters multiple primer sets at once.
#'
#' @param primers List with primer data frames.
#' @param templates List with template data frames.
#' @param active.constraints Strings giving the constraints that are to be checked.
#' @param settings List with settings.
#' @param updateProgress Progress callback function for shiny.
#' @return A list with filtered primer data frames.
#' @keywords internal
filter.comparison.primers <- function(primers, templates, active.constraints, 
    settings, updateProgress = NULL) {

    names <- names(primers)
    no.structures <- TRUE
    for (i in seq_along(primers)) {
        primer.df <- primers[[i]]
        template.df <- templates[[i]]
        mode.directionality <- get.analysis.mode(primer.df)
        filtered.data <- cascaded.filter.quick(primer.df, template.df, settings, 
            active.constraints, mode.directionality,
            active.constraints, no.structures, updateProgress)
        filtered.result <- filtered.data$data
        template.df <- update_template_cvg(template.df, primer.df, mode.directionality)
        primers[[i]] <- filtered.result
        templates[[i]] <- template.df
    }
    names(primers) <- names
    result <- list(primers = primers, templates = templates)
    #print("out len:")
    #print(length(result$primers))
    #print("out len templates:")
    #print(length(result$templates))
    return(result)
}
#' Filter a Set of Primers.
#'
#' Filters a primer set according to the constraints specified via 
#' \code{settings} and \code{active.constraints} such that all primers
#' that do not fulfill the constraints are removed from \code{primer.df}.
#'
#' @param primer.df A \code{Primers} object containing the primers
#' to be filtered.
#' @param template.df A \code{Templates} object with the template sequences 
#' that are to be covered by \code{primer.df}.
#' @param settings A \code{DesignSettings} object specifying the parameters 
#' for filtering the primers.
#' @param active.constraints The constraints that are to be used for filtering
#' primers. By default, \code{active.constraints} is set to \code{NULL} such that
#' all active constraints are used. 
#' @return A \code{Primers} object containing only those primers 
#' fulfilling all specified constraints.
#' @note Please note that some constraints can only be computed if additional software is installed,
#' please see \code{\link{DesignSettings}} for an overview.
#' @family primer functions
#' @export
#' @keywords Primers
#' @examples
#' data(Ippolito)
#' filename <- system.file("extdata", "settings", 
#'              "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
#' settings <- read_settings(filename)
#' # Only retain the primers fulfilling the GC clamp constraint:
#' filtered.df <- filter_primers(primer.df, template.df, settings,
#'                  active.constraints = c("gc_ratio"))
filter_primers <- function(primer.df, template.df, settings,
                    active.constraints = names(constraints(settings))) {
   
    # nb: only input_filtering_constraint settings are considered
    if (length(settings) == 0 || !is(settings, "DesignSettings")) {
        stop("Please provide a DesignSettings object.")
    }
    if (!is(primer.df, "Primers")) {
        stop("Please input a primer data frame.")
    }
    if (!is(template.df, "Templates")) {
        stop("Please provide a valid template data frame.")
    }
    filter.result <- cascaded.filter.quick(primer.df, template.df, settings, 
        to.compute.constraints = active.constraints, 
        mode.directionality = get.analysis.mode(primer.df), 
        active.constraints = active.constraints, 
        updateProgress = NULL)
    return(filter.result$data) # return remaining primers after filtering
}
#' Filter By Constraints
#' 
#' Remove primers that do not fulfill the current constraints (evaluate all primers).
#'
#' @param filtered.df Primer data frame.
#' @param constraint.df Data frame with constraint values.
#' @param current.constraints List with constraint settings.
#' @param active.constraints Strings giving the names of active constraints.
#' @param mode.directionality Direction of primers
#' @param template.df Template data frame.
#' @return A list containing the filtered primer data frame, as well as a 
#' data frame of the excluded primers and the used filtering settings.
#' @keywords internal
filter.by.constraints <- function(filtered.df, constraint.df, current.constraints, 
    active.constraints, mode.directionality = c("fw", "rev", "both"), template.df) {
    
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    #print(colnames(constraint.df))
    # evaluate currently active constraints
    eval.df <- eval.constraints(constraint.df, current.constraints, active.constraints, 
        mode.directionality, filtered.df)
    eval.t <- evaluate.template.constraints(constraint.df, current.constraints,
                                active.constraints, mode.directionality)
    # add hexamer coverage constraint if present
    eval.d <- lapply(eval.t, function(x) unlist(lapply(x, function(y) paste(y, collapse = ","))))
    eval.d <- do.call(cbind, eval.d)
    final.res <- merge.template.decisions(eval.t)
    idx.add <- which(!colnames(eval.df) %in% constraint.df)
    idx.rep <- which(colnames(eval.df) %in% constraint.df)
    add.df <- constraint.df
    if (length(idx.add) != 0) {
        add.df[, colnames(eval.df)[idx.add]] <- eval.df[, idx.add]
    }
    if (length(idx.rep) != 0) {
        add.df[, colnames(eval.df)[idx.rep]] <- eval.df[, idx.rep]
    }
    m <- match(colnames(add.df), colnames(filtered.df))
    # filtered.df becomes matrix here ...
    if (any(!is.na(m))) {
        filtered.df[, m[which(!is.na(m))]] <- add.df[, which(!is.na(m)), drop = FALSE]
    }
    if (any(is.na(m))) {
        filtered.df <- cbind(filtered.df, add.df[, which(is.na(m)), drop = FALSE])
    }
    # consider template-based constraints: remove covered seqs
    if ("hexamer_coverage" %in% active.constraints) {
        # don't filter coverage events with complementary 3' ends
        # adjust selection of coverage events: free-pass for hexa cvg
        hexa <- lapply(strsplit(constraint.df$hexamer_coverage, split = ","), function(x) as.logical(x))
        final.res <- lapply(seq_along(final.res), function(x) hexa[[x]] | final.res[[x]])
    }
    if (length(final.res) != 0) {
        # only update template-specific coverage data if we compute template-coverage specific constraints
        # determine indices of cvg events to be retained:
        sel <- lapply(final.res, function(x) which(x)) # selection of coverage events to retain according to cvg constraints
        filtered.df <- update.cvg.data(filtered.df, sel, template.df, mode = "on_target", active.constraints)
        if ("primer_coverage" %in% colnames(constraint.df)) {
            if ("constraints_passed_off_T" %in% colnames(constraint.df)) {
                # filtering of off-target binding events to update specificity
                off.sel <- lapply(strsplit(constraint.df$constraints_passed_off_T, split = ","), function(x) which(as.logical(x)))
            } else {
                # no off-target events filtered -> select all
                off.sel <- lapply(strsplit(constraint.df$Off_Covered_Seqs, split = ","), function(x) seq_len(length(x)))
            }
            filtered.df <- update.cvg.data(filtered.df, off.sel, template.df, mode = "off_target", active.constraints)
        }
    }
    # filter primers by considering all computed constraints
    passed <- rep(TRUE, nrow(filtered.df))
    if (length(eval.df) != 0) {
        passed <- apply(eval.df, 1, all)
    }
    # filter primer.df
    excluded <- NULL
    passed.idx <- which(passed)
    if (length(which(!passed) != 0)) {
        # primers were excluded
        excluded <- filtered.df[which(!passed), , drop = FALSE]
    }
    if (length(passed.idx) != 0) {
        filtered.df <- filtered.df[which(passed), ]
    } else {
        # nothing passed the constraints
        filtered.df <- filtered.df[FALSE, ]
    }
    result <- list(Filtered = filtered.df, Excluded = excluded, Passed_Idx = passed.idx)
    return(result)
}
#' Cascaded Filter
#'
#' Filter primers in a cascaded fashion. 
#'
#' At each constraint evaluation all primers that do not fulfill the current constraint are removed.
#'
#' Constraints that are specified in \code{to.compute.constraints} are computed on the fly.
#'
#' @param primer.df Primer data frame.
#' @param template.df Template data frame.
#' @param settings Settings object.
#' @param to.compute.constraints Names of constraints that still have to be computed.
#' @param mode.directionality Primer direction.
#' @param active.constraints The constraints that are to be used for filtering.
#' If \code{active.constraints} is \code{NULL}, all filtering constraints are used.
#' @param no.structures Whether dimerization structures shall be computed.
#' @param updateProgress Progress callback function for shiny.
#' @return The filtered primer data frame.
#' @keywords internal
cascaded.filter.quick <- function(primer.df, template.df, settings, 
    to.compute.constraints, 
    mode.directionality = c("fw", "rev", "both"), 
    active.constraints = NULL, no.structures = FALSE, updateProgress = NULL) {

    if (length(primer.df) == 0 || length(template.df) == 0) {
        return(NULL)
    }
    if (!is(primer.df, "Primers")) {
        stop("Please supply a 'Primers' object.")
    }
    if (!is(template.df, "Templates")) {
        stop("Please supply a 'Templates' object.")
    }
    if (!is(settings, "DesignSettings")) {
        stop("Please supply a 'DesignSettings' object.")
    }
    PCR_conditions <- PCR(settings)
    other.settings <- conOptions(settings)
    allowed.region.definition <- other.settings$allowed_region_definition
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    filtered.df <- primer.df
    if (is.null(active.constraints)) {
        # no constraints provided -> use all constraints to compute
        active.constraints <- to.compute.constraints
    }
    # use specified constraints 
    current.constraints <- constraints(settings)
    m <- match(active.constraints, names(current.constraints))
    if (any(is.na(m))) {
        # couldn't find an input constraint
        stop(paste("Unknown 'active.constraints' specified.",
            "The constraints should be an available constraint from the input",
            "settings object. Unknown constraints are: ",
            paste(active.constraints[is.na(m)], collapse = ",")))
    }
    current.constraints <- current.constraints[m]
    if (length(current.constraints) == 0) {
        return(list(data = primer.df, stats = NULL, excluded = NULL, used_settings = settings))  # nothing to filter 
    }
    stat.N <- length(current.constraints)  # +1 for initial nbr of primers
    stat.df <- data.frame(Constraint = names(current.constraints), Direction = mode.directionality, Remaining = rep(0, stat.N), Excluded = rep(0, 
        stat.N), Current_Coverage = rep(0, 
        stat.N), Time = rep(0, stat.N), stringsAsFactors = FALSE)  # determine statistics about filtering
    max.cvg <- get_cvg_ratio(filtered.df, template.df)
    # load settings from input data
    excluded.df <- primer.df[FALSE, ]
    last.constraints <- current.constraints
    for (i in seq_along(current.constraints)) {
        active.constraints <- names(current.constraints)[i]
        #print(active.constraints)
        ####### Time measurement for filtering step
        ptm <- proc.time()
        if (active.constraints %in% to.compute.constraints) {
            # value should be computed
            constraint.df <- compute.constraints(filtered.df, mode.directionality, 
                template.df, settings, active.constraints, 
                no.structures = no.structures,
                updateProgress = updateProgress)
        } else {
            constraint.df <- filtered.df  # value isn't recomputed -> should be in filtered.df
        }
        ptm <- proc.time() - ptm
        ####### Time measurement for filtering step ends
        filtered.df.pure <- filtered.df  # wasn't filtered yet :-)
        # evaluate currently active constraints
        filter.result <- filter.by.constraints(filtered.df, constraint.df, current.constraints, 
            active.constraints, mode.directionality, template.df)
        filtered.df <- filter.result$Filtered
        excluded <- filter.result$Excluded
        if (length(excluded) != 0) {
            excluded.df <- my_rbind(excluded.df, excluded)
        }
        max.cvg <- get_cvg_ratio(filtered.df, template.df)
        stat.df[i, "Constraint"] <- active.constraints
        stat.df[i, "Current_Coverage"] <- max.cvg
        nbr.excluded <- ifelse(length(excluded) == 0, 0, nrow(excluded))
        stat.df[i, "Excluded"] <- nbr.excluded
        stat.df[i, "Remaining"] <- nrow(filtered.df)  # included/remaining observations
        stat.df[i, "Time"] <- ptm["elapsed"]
        percentage.remaining <- paste0(round(nrow(filtered.df) / nrow(primer.df) * 100, 2), "%")
        message("\to ", active.constraints, ": ", nrow(filtered.df), " primers remaining (", percentage.remaining, ")")
        if (nrow(filtered.df) == 0) {
            # nothing more to filter.
            break
        }
    }
    s <- 1  # start with first entry now (no more initial)
    e <- length(stat.df$Constraint)
    excluded.df$Exclusion_Reason <- unlist(lapply(s:e, function(x) rep(stat.df$Constraint[x], 
        stat.df$Excluded[x])))
    # does excluded.df have some colnames missing?
    diff <- setdiff(colnames(filtered.df), colnames(excluded.df))
    if (length(diff) != 0 && nrow(excluded.df) != 0) {
        excluded.df[, diff] <- rep(NA, nrow(excluded.df))
    }
    # update melting_temp_diff
    if (length(filtered.df) != 0 && "melting_temp_diff" %in% colnames(filtered.df) && nrow(filtered.df) != 0) {
        temp.diff <- get.melting.temp.diff(filtered.df$Tm_C_fw, filtered.df$Tm_C_rev)
        filtered.df[, "melting_temp_diff"] <- temp.diff
    }
    stat.df$Excluded_Percentage <- stat.df$Excluded/(stat.df$Remaining + stat.df$Excluded)
    result <- list(data = filtered.df, stats = stat.df, excluded = excluded.df)
    return(result)
}

#' Writes Filtering Data Sets to Disk.
#' @param filtered.df A filtered \code{Primers} set.
#' @param excluded.df A set of \code{Primers} that were excluded.
#' @param results.loc The location where to store the data.
#' @param tag A tag for the output files.
#' @param stat.df Data frame with statistics of the filtering procedure.
#' @param settings A \code{DesignSettings} object.
#' @return No return value, writes output to disk.
#' @keywords internal
store.filtering.sets <- function(filtered.df, excluded.df, results.loc, tag = "", stat.df = NULL, settings = NULL) {
    if (length(results.loc) != 0) {
        write.csv(excluded.df, file = file.path(results.loc, paste0("excluded_primers_", tag, ".csv")), 
            row.names = FALSE)
        write.csv(filtered.df, file = file.path(results.loc, paste0("filtered_primers_", tag, ".csv")), 
            row.names = FALSE)
        if (length(stat.df) != 0) {
            write.csv(stat.df, file = file.path(results.loc, paste0("filtered_primers_stats_", tag, ".csv")), 
                row.names = FALSE)
        }
        if (length(settings) != 0) {
            write_settings(settings, file.path(results.loc, paste0("relaxed_settings_", tag, ".xml")))
        }
    }
    invisible()
}
#' Filtering for the Optimization
#' 
#' Filter primers according to constraints and relax constraints if necessary.
#'
#' Constraints are relaxed if the \code{required.cvg} could not be reached with the input constraints.
#' 
#' @param primer.df Primer data frame.
#' @param template.df Template data frame.
#' @param settings Settings object.
#' @param mode.directionality Primer direction.
#' @param required.cvg Required ratio of covered templates.
#' If \code{required.cvg} is set to 0, the constraints are not relaxed.
#' @param target.temps Target melting temperature of the primers in Celsius.
#' This argument is only required if we try to match the melting temperatures of another primer set,
#' e.g. when first optimizing forward and then optimizing reverse primers.
#' @param updateProgress Progress callback function for shiny.
#' @param results.loc Directory where the filtering results should be stored.
#' @return The filtered primer data frame with respect to \code{required.cvg}.
#' @keywords internal
cascaded.filter <- function(primer.df, template.df, settings, mode.directionality = c("fw", "rev", "both"), 
                            required.cvg = 1, target.temps = NULL, 
                            updateProgress = NULL, results.loc = NULL) {
    
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    constraint.settings <- filters(settings)
    to.compute.constraints <- names(constraint.settings)  # compute all constraints given as input
    message("a) Filtering of primers")
    # set no.structures to TRUE since this is the design algorithm
    filtered.data <- cascaded.filter.quick(primer.df, template.df, settings, 
        to.compute.constraints, mode.directionality, no.structures = TRUE,
        updateProgress = updateProgress)  # filter without relaxation
    if (length(filtered.data$excluded) == 0) {
        return(list(excluded = filtered.data$excluded, data = filtered.data$data, 
            stats = filtered.data$stats, used_settings = settings, 
            stats_relaxed = NULL))
    }
    excluded.df <- filtered.data$excluded
    filtered.df <- filtered.data$data
    stat.df <- filtered.data$stats
    # relax here:
    stat.df$Excluded_Percentage <- stat.df$Excluded/(stat.df$Remaining + stat.df$Excluded)
    # save results up to now
    store.filtering.sets(filtered.df, excluded.df, results.loc, tag = "pre", stat.df = stat.df)
    message("b) Relaxation of constraints")
    # relax constraints:
    relaxed.data <- relax.constraints(settings, 
        filtered.df, excluded.df, stat.df, template.df, mode.directionality, 
        required.cvg, target.temps = target.temps,
        results.loc = results.loc)
    excluded.df <- relaxed.data$excluded
    filtered.df <- relaxed.data$filtered  # included seqs
    stats.relax <- relaxed.data$stats
    relaxed.settings <- relaxed.data$new_settings
    #max.cvg <- get_cvg_ratio(filtered.df, template.df)
    result <- list(data = filtered.df, stats = stat.df, excluded = excluded.df, 
                   used_settings = relaxed.settings, stats_relax = stats.relax)
    return(result)
}
#' Uncovered Templates.
#'
#' Computes a data frame containing the templates that are
#' not yet covered.
#'
#' @param filtered.df An object of class \code{Primers}.
#' @param template.df An object of class \code{Templates}.
#' @return A Templates data frame containing the missing templates.
#' @keywords internal
get.missing.df <- function(filtered.df, template.df, Tm.brackets, settings, mode.directionality) {
    # TODO: missing should be based on the temperature primer sets and not all primers together ..!! TODO
    # build df of templates we still haven't covered
    if (!"Covered_Seqs" %in% colnames(filtered.df)) {
        covered.idx <- NULL
    } else {
        # get covered templates per temperature-set
        cvg.info <- get_max_set_coverage(filtered.df, template.df, Tm.brackets, settings, mode.directionality, max.only = FALSE)
        covered.idx <- lapply(cvg.info, function(x) match(attr(x,"covered_templates"), template.df$ID))
    }
    # determine primer candidates that can still improve the coverage
    missing.idx <- NULL
    if (length(covered.idx) == 0) {
        missing.idx <- 1:nrow(template.df)
    } else {
        to.cover <- 1:nrow(template.df)
        # consider all uncovered templates from all temperature sets:
        missing.idx <- unique(unlist(lapply(covered.idx, function(x) setdiff(to.cover, x))))
    }
    missing.df <- template.df[missing.idx, ]
    return(missing.df)
}
#' Computation of Coverage Gain.
#'
#' Computes the coverage gain from \code{covered.seqs}.
#'
#' @param covered.seqs List with covered sequences.
#' @param template.df A \code{Templates} data frame.
#' @param missing.df A \code{Templates} data frame 
#' containing only the templates that still need to be covered.
#' @param candidate.df A \code{Primers} data frame containing
#' candidate primers.
#' @param con.names The constraint to evaluate the coverage gain for
#' upon being relaxed.
#' @param constraint.limits A list with constraint limits.
#' @param feasible.only Whether only feasible coverage gains 
#' are to be outtputed. Here, \emph{feasible} relates to 
#' coverage gains that can be obtained directly with the next relaxation.
#' @return The number of covered sequences to be gained.
#' @keywords internal
get.cvg.gain <- function(covered.seqs, template.df, missing.df, candidate.df, 
                        con.names, constraint.limits, feasible.only = FALSE) {
    # identify the coverage for covered.seqs for a constraint
    # covered seqs: identifiers of covered seqs
    con.name <- unique(con.names)
    cur.limit <- constraint.limits[[con.name]]
    values <- get.constraint.values(unique(con.names), candidate.df, "both")
    fw.idx <- which(candidate.df$Forward != "")
    rev.idx <- which(candidate.df$Reverse != "")
    covered.idx <- covered.seqs.to.idx(covered.seqs, template.df)
    if (feasible.only) {
        feasible <- apply.constraint(values$fw, values$rev, cur.limit["min"], 
                         cur.limit["max"], fw.idx, rev.idx, "both")
        sel <- which(apply(feasible, 1, all))
        covered.idx <- covered.idx[sel]
    }
    unique.covered.idx <- unique(unlist(lapply(covered.idx, function(x) intersect(missing.df$Identifier, template.df[x, "Identifier"]))))
    cvg <- length(unique.covered.idx)
    return(cvg)
}
#' Gain of Coverage by Excluded Primers.
#'
#' Computes a data frame on the excluded sequences per constraint.
#'
#' @param candidate.df An object of class \code{Primers} containing
#' excluded primers that are considered for addition to \code{filtered.df}.
#' @param constraint.settings A list with the current constraint settings.
#' @param constraint.limits The current constraint limits.
#' @param relax.df Data frame with count of relaxations per constraint.
#' @return A data frame with exclusion data.
#' @keywords internal
build.gain.df <- function(candidate.df, constraint.settings, constraint.limits, relax.df) {
    if (length(candidate.df) == 0 || nrow(candidate.df) == 0) {
        # nothing to be added
        return(NULL)
    }
    # determine cvg gain of primer candidates per exclusion reason
    df <- relax.df
    if (nrow(df) != 0) {
        # only consider constraints for relaxation that are associated with an exclusion reason
        idx <- unlist(lapply(df$Constraint, function(x) length(which(candidate.df$Exclusion_Reason == x))))
        sel.idx <- which(idx != 0)
        # select 'unrelaxed' constraints first and then select by recommended order of relaxation (relax unimportant constraints first)
        relaxable.idx <- which(unlist(lapply(df$Constraint, function(x) 
                                   any(constraint.settings[[x]] != constraint.limits[[x]]))))
        idx <- intersect(sel.idx, relaxable.idx)
        df <- df[idx,]
        df <- df[order(df$Relax_Count, match(df$Constraint, getOption("openPrimeR.relax_order"))), ]
    }
    #message("Coverage gain:")
    #message(df)
    return(df)
}
#' Determination of Maximal Coverage.
#'
#' Determines the maximal coverage ratio of a set of primers 
#' for primer subsets valid for a certain temperature range.
#' a certain melting temperature range.
#'
#' @param primer.df An object of class \code{Primers}.
#' @param template.df An objectc of class \code{Templates}.
#' @param Tm.brackets A data frame with temperature information.
#' @param settings A \code{DesignSettings} object.
#' @param mode.directionality The direction of the primers.
#' @param max.only Whether only the maxium coverage shall be returned.
#' If \code{max.only} is \code{FALSE}, the coverage ratios of all 
#' melting temperature sets according to \code{Tm.brackets} are returned.
#' @return The maximal coverage ratio of a primer set if \code{max.only} is \code{TRUE}
#' or the coverages of all melting temperature seets if \code{max.only} is \code{FALSE}.
#' @keywords internal
get_max_set_coverage <- function(primer.df, template.df, Tm.brackets, settings, mode.directionality, max.only = TRUE) {
    # determine coverage for individual melting temperature ranges to ensure that we have a set with max cvg that is temperature compatible!
    #print("get max set coverage")
    PCR.settings <- PCR(settings)
    other.settings <- conOptions(settings)
    allowed.region.definition <- other.settings$allowed_region_definition
    Tm.sets <- compute.Tm.sets(primer.df, template.df, Tm.brackets, settings, 
        mode.directionality, PCR.settings$primer_concentration, PCR.settings$template_concentration,
        PCR.settings$Na_concentration, PCR.settings$Mg_concentration, 
        PCR.settings$K_concentration, PCR.settings$Tris_concentration,
        other.settings$allowed_mismatches, 
        other.settings$allowed_other_binding_ratio, other.settings$allowed_stop_codons, allowed.region.definition, 
        other.settings$disallowed_mismatch_pos, FALSE, NULL, NULL, NULL, NULL)$sets  # get temperature-filtered sets
    cvg.vals <- lapply(Tm.sets, function(x) get_cvg_ratio(x, template.df)) # get cvg of every Tm set
    if (max.only) {
        cvg.vals <- max(unlist(cvg.vals))
    }
    #print("max set coverage done")
    return(cvg.vals)
}
#' Relaxation of Constraint Limits.
#'
#' Relaxes the constraint limits by moving according to the difference
#' between the \code{initial.limits} and \code{initial.constraints}.
#'
#' @param cur.limits List with current constraint settings.
#' @param initial.limits List with initial coonstraint limits.
#' @param initial.constraints List with initial constraint settings before relaxing.
#' @param con.name The constraint for which the settings are to be changed.
#' @return A list with relaxed constraint limits.
#' @keywords internal
set.new.limits <- function(cur.limits, initial.limits, initial.constraints, con.name = NULL) {
    # relaxes the constraint limits
    # cur.limit: current constraint limits
    # initial.limits: limits for constraints initially
    # initial.constraints: settings for constraints initially

    # 1. Determine the differences between the initial.constraints and the initial.limits 
    # -> determines the speed of the relaxation.
    if (length(con.name) == 0) {
        con.name <- names(cur.limits)
    }
    constraint.delta <- lapply(seq_along(con.name), function(x) initial.limits[[con.name[x]]] - initial.constraints[[con.name[x]]])
    names(constraint.delta) <- con.name
    # 2. update all limits according to the shift from the initial constraints
    new.limits <- cur.limits
    m <- match(con.name, names(new.limits))
    new.limits[con.name] <- lapply(seq_along(con.name), function(x) new.limits[[con.name[x]]] + constraint.delta[[con.name[x]]])
    names(new.limits)[m] <- con.name
    return(new.limits)
}
#' Get the Values of a Constraint.
#'
#' @param con.name The name of the constraint.
#' @param cur.candidates The \code{Primers} data frame where the values should be retrieved.
#' @param mode.directionality The direction for which values should be retrieved.
#' @return The constraint values corresponding to \code{con.name} for the primers \code{cur.candidates}.
#' @keywords internal
get.constraint.values <- function(con.name, cur.candidates, mode.directionality) {
    value.idx <- unlist(get.constraint.value.idx(con.name, cur.candidates))  # idx of constraint values in constraint.df
    val.name <- colnames(cur.candidates)[value.idx]
    values <- asS3(cur.candidates)[, val.name, drop = FALSE]
    rev.col <- grep("_rev", colnames(values))  # index of rev col value
    # for both: retrieve fw & rev!
    if (mode.directionality == "both") {
        modes <- c("fw", "rev")
    } else {
        modes <- mode.directionality
    }
    result <- vector("list", 2)
    names(result) <- c("fw", "rev")
    for (i in seq_along(modes)) {
        mode <- modes[i]
        col.id <- grep(paste("_", mode, sep = ""), colnames(values))  # index of fw col value
        if (length(col.id) != 0) {
            # constraint with fw & rev values
            out.values <- as.numeric(values[, col.id])
        } else {
            # constraint with a single value for both directions
            out.values <- values[, 1]
        }
        result[[mode]] <- out.values
    }
    return(result)
}
#' Augmentation of Primer Coverage.
#'
#' Computes the coverage for the primers in \code{primer.df} 
#' that is still missing such that the relaxation procedure can adjust
#' appropriate constraints.
#'
#' @param primer.df A \code{Primers} object for which
#' the primer coverage shall be augmented.
#' @param template.df A \code{Templates} object.
#' @param settings A \code{DesignSettings} object giving
#' the parameters for coverage computations.
#' @param partial Whether all missing primer coverage values should be 
#' computed. If \code{partial} is \code{TRUE}, only 
#' the coverage values of the primers that were excluded due to
#' the specified \code{constraint} are computed.
#' @param constraint A character vector specifying the exclusion
#' reason for which the partial augmentation should take place.
#' @return A \code{Primers} object with augmented coverage entries.
#' @keywords internal
augment.primer.cvg <- function(primer.df, template.df, settings, 
                               partial = FALSE, constraint = NULL) {
    if (!is(primer.df, "Primers")) {
        stop("Please input a valid primer data frame.")
    }
    if (!is(template.df, "Templates")) {
        stop("Please input a valid template data frame.")
    }
    if (!is(settings, "DesignSettings")) {
        stop("Please input a valid settings object.")
    }
    is.exclusion.df <- "Exclusion_Reason" %in% colnames(primer.df)
    if (!"primer_coverage" %in% colnames(primer.df)) {
        primer.df$primer_coverage <- NA
        primer.df$Covered_Seqs <- NA
    }
    # compute only missing coverage entries
    idx <- which(is.na(primer.df$primer_coverage))
    ori.idx <- idx
    if (partial && is.exclusion.df && length(constraint) != 0) {
        con.idx <- which(primer.df$Exclusion_Reason == constraint)
        idx <- intersect(idx, con.idx) # primers excluded for 'constraint' with missing cvg
    }
    if (length(idx) == 0) {
        # nothing to adjust
        return(primer.df)
    }
    con.s <- ifelse(length(constraint) != 0, paste0(constraint, ": "), "")
    message("\to ", con.s, "Augmenting primer coverage for ", length(idx), " of ", length(ori.idx), " target primers.")
    augment.df <- check_constraints(primer.df[idx,], template.df, settings, "primer_coverage")
    existing.cols <- intersect(colnames(primer.df), colnames(augment.df))
    primer.df[idx,existing.cols] <- augment.df[, existing.cols]
    new.cols <- setdiff(colnames(augment.df), colnames(primer.df))
    for (i in seq_along(new.cols)) {
        col <- new.cols[i]
        primer.df[idx, col] <- augment.df[idx, col]
    }
    return(primer.df)
}
#' Relaxation of Constraints
#'
#' Relax constraints trying to reach the target coverage ratio.
#'
#' @param settings A \code{DesignSettings} object.
#' @param filtered.df Data set of primers that fulfilled all constraints of the filtering procedure.
#' @param excluded.df Data frame with excluded primers from the first filtering round.
#' @param stat.df Data frame with statistics from filtering.
#' @param template.df Template data frame.
#' @param mode.directionality Primer direction.
#' @param required.cvg Required ratio of covered templates.
#' @param target.temps Target melting temperature values.
#' @param results.loc The location where intermediary results should be stored.
#' @return A list containing information about the relaxation as well as the filtered primers.
#' @keywords internal
relax.constraints <- function(settings, filtered.df, excluded.df, stat.df, template.df, 
                              mode.directionality = c("fw", "rev"), required.cvg, 
                              target.temps = NULL,
                              results.loc = NULL) {

    if (length(results.loc) != 0) {
        # create result location
        results.loc <- file.path(results.loc, "relaxation")
        dir.create(results.loc, showWarnings = FALSE)
    }
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    initial.constraints <- filters(settings)
    active.constraints <- names(initial.constraints)
    initial.limits <- filterLimits(settings)
    initial.opti.constraints <- opti(settings)
    initial.opti.limits <- optiLimits(settings)
    new.filtered.df <- filtered.df  # updated filtered data from relaxation
    if (nrow(new.filtered.df) != 0) {
        new.filtered.df$Exclusion_Reason <- NA # add exclusion reason column
    }
    new.excluded.df <- excluded.df
    Tm.brackets <- create.Tm.brackets(new.filtered.df, template.df, settings, target.temps)
    new.filtered.df <- Tm.brackets$primers
    max.cvg <- get_max_set_coverage(new.filtered.df, template.df, Tm.brackets, settings, mode.directionality)
    if (is.na(max.cvg)) {
        stop("The maximal coverage is 'NA', which indicates that ",
             "primer coverage was not available after filtering.",
             "This shouldn't have had happened, sorry.")
    } 
    while.count <- 0
    relax.df <- data.frame("Constraint" = names(initial.constraints), "Relax_Count" = 0, stringsAsFactors=FALSE)
    # count the total gained coverage:
    total.gain.df <- data.frame(Exclusion_Reason = active.constraints,
                                Actual_Gain = 0, stringsAsFactors=FALSE)
    target.cvg <- required.cvg
    # Augment primer coverage for all primers (not only for missing templates since runtime would increase for multiple cvg calls, would need to recompute later)
    #t <- proc.time()["elapsed"]
    missing.df <- NULL
    if (max.cvg < target.cvg) {
        new.excluded.df <- augment.primer.cvg(excluded.df, template.df, settings)
        missing.df <- get.missing.df(new.filtered.df, template.df, Tm.brackets, settings, mode.directionality)
    }
    #t <- proc.time()["elapsed"] - t
    #print(t)
    while (max.cvg < target.cvg) {
        if (length(results.loc) != 0) {
            # store relaxation results
            store.filtering.sets(new.filtered.df, new.excluded.df, results.loc, tag = as.character(while.count), stat.df = relax.df, settings = settings)
        }
        while.count <- while.count + 1
        # improve coverage as long as we haven't reached the target
        target.cvg <- required.cvg # TODO: change to feasible coverage if we have all the coverage values computed...
        if ("Covered_Seqs" %in% colnames(new.filtered.df)) {
            prev.cvg <- unique(unlist(strsplit(new.filtered.df$Covered_Seqs, split = ",")))
        } else {
            # everything was filtered before cvg was computed.
            prev.cvg <- NULL
        }
        # 1. Build Exclusion Reason data frame
        #message(paste("Relaxation run: ", while.count, sep = ""))
        gain.df <- build.gain.df(new.excluded.df, filters(settings), filterLimits(settings), relax.df)
        ######################
        print("#########")
        print(gain.df)
        print("#########")
        #print("Melting temperatures of excluded seqs: ")
        #print(new.excluded.df$melting_temp)
        ########################
        if (length(gain.df) == 0 || nrow(gain.df) == 0) {
            # there's no constraint to relax anymore to bring forth an increase in coverage
            # DEBUG:
            #cvg.to.gain <- unique(unlist(strsplit(new.excluded.df$Covered_Seqs, split = ",")))
            #cur.cvg <- unique(unlist(strsplit(new.filtered.df$Covered_Seqs, split = ",")))
            #possible.gain <- setdiff(cvg.to.gain, cur.cvg)
            if ("melting_temp_diff" %in% names(opti(settings)) && "melting_temp" %in% colnames(new.filtered.df) && 
                any(constraints(settings)$melting_temp_diff != constraintLimits(settings)$melting_temp_diff)) {
                    max.steps <- 10 # limit the number of steps, in case it's just not possible (e.g. cannot relax anymore due to user input)
                    step <- 0
                    # this procedure uses only 'new.filtered.df', assuming that we couldn't relax anymore and it's only a matter of 'Tm diff' to reach the target.
                    while (max.cvg < target.cvg && step < max.steps) {
                        step <- step + 1
                        message("Final melting temp diff relaxation ...")
                        # adjust melting temp diff until we reach target cvg
                        constraints(settings)$melting_temp_diff <- relax.opti.constraints(opti(settings), 
                                                                   initial.opti.limits, initial.opti.constraints)$melting_temp_diff
                        constraintLimits(settings)$melting_temp_diff <- set.new.limits(optiLimits(settings), initial.opti.limits, 
                                                                    initial.opti.constraints, "melting_temp_diff")$melting_temp_diff
                        Tm.brackets <- create.Tm.brackets(new.filtered.df, template.df, settings, target.temps)
                        new.filtered.df <- Tm.brackets$primers
                        max.cvg <- get_max_set_coverage(new.filtered.df, template.df, Tm.brackets, settings, mode.directionality)
                    }
            }
            break
        } 
        # select the first constraint in gain.df: most promising exclusion reason for relaxing
        con.name <- gain.df$Constraint[1] # gain.df is sorted in the order we want to relax constraints in
        message("Relaxing constraint: ", con.name)
        con.idx <- which(active.constraints == con.name)
        # determine candidate primers 
        # only re-check constraints for those primers that give us a gain in coverage
        cvd.idx <- covered.seqs.to.idx(new.excluded.df$Covered_Seqs, template.df)
        cvg.idx <- which(sapply(cvd.idx, function(x) length(intersect(missing.df$Identifier, template.df$Identifier[x])) != 0))
        new.excluded.df <- new.excluded.df[cvg.idx, ]
        exclusion.idx <- which(new.excluded.df$Exclusion_Reason == con.name)
        #message("Candidate idx: ", paste(idx, collapse = ","))
        # augment primer coverage for all primers that failed the current constraint
        cur.candidates <- new.excluded.df[exclusion.idx,]
        # determine constraints that still need to be computed and checked
        to.compute.constraints <- active.constraints[tail(seq_along(filters(settings)), 
            length(filters(settings)) - con.idx)]  # constraints that need to be computed (all later constraints than con.idx)
        active.constraints.check <- active.constraints[tail(seq_along(filters(settings)), 
            length(filters(settings)) - con.idx + 1)]  # constraints to be checked for fulfillment (all later constraints including the current constraint)
        # never compute cvg!
        to.compute.constraints <- to.compute.constraints[!to.compute.constraints %in% c("primer_coverage", "primer_specificity")]
        message("Computing constraints: ", paste(to.compute.constraints, collapse = ","))
        # set new constraint value to limits
        new.setting <- constraintLimits(settings)[[con.name]]
        # relax the constraint
        message("\to ", con.name, ": Relaxed constraint from [", 
            paste0(round(filters(settings)[[con.name]], 2), 
                    collapse = ", "), 
             "] to [", 
              paste0(round(new.setting, 2), collapse = ","), 
              "].")
        relax.df[relax.df$Constraint == con.name, "Relax_Count"] <- relax.df[relax.df$Constraint == con.name, "Relax_Count"] + 1
        constraints(settings)[[con.name]] <- new.setting  # update constraints
        # evaluate constraints again
        # set no.structures to TRUE since we are designing
        filter.result <- suppressMessages(cascaded.filter.quick(cur.candidates, template.df, settings, 
            to.compute.constraints = to.compute.constraints, 
            mode.directionality = mode.directionality,
            active.constraints = active.constraints.check, 
            no.structures = TRUE, updateProgress = NULL))
        # update the filtered primers
        new.primers <- filter.result$data
        new.filtered.df <- my_rbind(new.filtered.df, new.primers)
        # update excluded primers:
        if (length(exclusion.idx) != 0) {
            # remove old entries
            new.excluded.df <- new.excluded.df[-exclusion.idx, ]
            new.excluded.df <- my_rbind(new.excluded.df, filter.result$excluded)
        }
        #print(new.excluded.df) # many things are NA here, why?
        # store actual coverage gain
        if (con.name %in% new.filtered.df$Exclusion_Reason) {
            new.covered <- unique(unlist(strsplit(new.filtered.df$Covered_Seqs[new.filtered.df$Exclusion_Reason == con.name], split = ",")))
        } else {
            new.covered <- NULL
        }
        additional.cvg <- length(setdiff(new.covered, prev.cvg))
        # store feasible and actual cvg gains
        total.gain.entry <- total.gain.df[total.gain.df$Exclusion_Reason == con.name, "Actual_Gain"]
        total.gain.entry <- total.gain.entry + additional.cvg
        total.gain.df[total.gain.df$Exclusion_Reason == con.name, "Actual_Gain"] <- total.gain.entry
        # update tm brackets to make sure that new melting temp data is taken into account
        Tm.brackets <- create.Tm.brackets(new.filtered.df, template.df, settings, target.temps)
        new.filtered.df <- Tm.brackets$primers
        # update cvg for while loop
        max.cvg <- get_max_set_coverage(new.filtered.df, template.df, Tm.brackets, settings, mode.directionality)
        message(paste("Current coverage is: ", max.cvg, sep = ""))
        # modify the limit of the relaxed constraint for possibly further relaxations
        constraintLimits(settings)[con.name] <- set.new.limits(filterLimits(settings), initial.limits, initial.constraints, con.name)[con.name]
        # update missing templates
        missing.df <- get.missing.df(new.filtered.df, template.df, Tm.brackets, settings, mode.directionality)
    }
    if ("melting_temp_diff" %in% names(constraints(settings))) {
        # keep only temperature primer sets with cvg >= target cvg
        Tm.cvg.vals <- unlist(get_max_set_coverage(new.filtered.df, template.df, Tm.brackets, settings, mode.directionality, max.only = FALSE))
        set.keep.idx <- which(Tm.cvg.vals >= target.cvg)
        # get melting temperature range of primers to keep
        Tm.df <- Tm.brackets$df
        sel.df <- Tm.df[set.keep.idx, c("min_Tm", "max_Tm")]
        min.Tm <- min(sel.df[,1])
        max.Tm <- max(sel.df[,2])
        # only retain primers between min.Tm and max.Tm, the others are not part of a set with neccesary cvg
        excluded.idx <- which(new.filtered.df$melting_temp < min.Tm | new.filtered.df$melting_temp > max.Tm)
        if (length(excluded.idx) != 0) {
            Tm.exclusion <- filtered.df[excluded.idx,]
            Tm.exclusion$Exclusion_Reason <- "melting_temp_range"
            new.excluded.df <- my_rbind(new.excluded.df, Tm.exclusion)
            new.filtered.df <- new.filtered.df[-excluded.idx,]
            # update current cvg (should be the same as before ...)
            max.cvg <- get_max_set_coverage(new.filtered.df, template.df, Tm.brackets, settings, mode.directionality)
        }
    }
    # store final set:
    store.filtering.sets(new.filtered.df, new.excluded.df, results.loc, tag = as.character(while.count), stat.df = relax.df, settings = settings)
    # store new statistics:
    excl.counts <- plyr::ddply(excluded.df, "Exclusion_Reason", plyr::summarize, "Count" = length(substitute(Exclusion_Reason)))
    # for additional cvg: only the last failed exclusion reason is considered
    counts <- plyr::ddply(new.filtered.df, "Exclusion_Reason", plyr::summarize, "Count" = length(substitute(Exclusion_Reason)), "Covered" = paste(unique(unlist(strsplit(substitute(Covered_Seqs), split = ","))), collapse = ","))
    if (nrow(counts) != 0) {
        # remove count entries that didn't fail a constraint
        na.idx <- which(is.na(counts$Exclusion_Reason ))
        if (length(na.idx) != 0) {
            counts <- counts[-na.idx,]
        }
    }
    new.stat.df <- data.frame(Constraint = active.constraints,
            Direction = mode.directionality, Remaining = 0,
            Time = 0, Excluded = 0, # nb time not stored here
            Current_Coverage = NA, 
            Excluded_Percentage = NA,
            stringsAsFactors = FALSE)
    m.c <- match(counts$Exclusion_Reason, new.stat.df$Constraint)
    new.stat.df$Remaining[m.c] <- counts$Count # Additional inclusions per category
    m.e <- match(excl.counts$Exclusion_Reason, new.stat.df$Constraint)
    new.stat.df$Excluded[m.e] <- excl.counts$Count # Previous exclusion counts
    new.stat.df$Excluded_Percentage <- new.stat.df$Excluded / (new.stat.df$Remaining + new.stat.df$Excluded) 
    # store current coverage: gain of cvg by every constraint
    new.stat.df$Current_Coverage <- total.gain.df$Actual_Gain / nrow(template.df)
    # return result:
    cvg.s <- paste0(round(max.cvg *100, 2), "%")
    message("-> Available number of candidates after relaxation: ", nrow(new.filtered.df), " (", cvg.s, " coverage)")
    out <- list(new_settings = settings, filtered = new.filtered.df, 
        excluded = new.excluded.df, stats = new.stat.df)
    return(out)
}

#' Update Constraint Settings.
#'
#' Updates the constraint settings with new values. Sets to the maximal
#' observed values within the limits or less if \code{relax.speed} is less than 1.
#'
#' @param values Observed constraint values considered for the update.
#' @param relax.speed The speed at which the constraints should be relaxed.
#' This value should be in the interval [0,1].
#' @param cur.limits List with current relaxation limits.
#' @param cur.setting List with current constraint settings.
#' @return Relaxed constraint settings according to the given \code{values} and \code{cur.limits}.
#' @keywords internal
set.new.constraint.value <- function(values, relax.speed, cur.limits, cur.setting) {
     #values <- get.constraint.values(con.name, cur.candidates, mode.directionality)[[mode.directionality]]

    # remove NA's from values
    values <- values[!is.na(values)]
    if (length(values) == 0) {
        return(NULL)        
    }
    new.setting <- cur.setting
    if (length(cur.limits) == 2) {
        sel.idx <- which(values >= cur.limits["min"] & values <= cur.limits["max"])  # only consider values that can be relaxed
        val <- values[sel.idx]
        if (length(val) != 0) {
            values.low <- val[val <= cur.setting["min"]]
            values.high <- val[val >= cur.setting["max"]]
            if (length(values.low) != 0) {
                new.setting["min"] <- unname(quantile(values.low, 1 - relax.speed))
            }
            if (length(values.high) != 0) {
                new.setting["max"] <- unname(quantile(values.high, relax.speed))
            }
            if (length(values.low) == 0 && length(values.high) == 0) {
                # it's crucial that even though there's nothing relaxable that the function returns something in this case for the relaxation procedure to work.
                return(cur.limits)
            }
            # don't set to limits when nothing changes but rather set to min/max of remaining
            # values to keep the constraints tight
            if (new.setting["min"] == cur.setting["min"] && new.setting["max"] == cur.setting["max"]) {
                if (length(values.low) != 0) {
                  # only change when there are still values
                  new.setting["min"] <- min(values.low, na.rm = TRUE)
                }
                if (length(values.high) != 0) {
                  new.setting["max"] <- max(values.high, na.rm = TRUE)
                }
            }
        } else {
            # no values within the current limits -> set constraints to the limits
            return(cur.limits)
        }
    } else if (length(cur.limits) == 1) {
        if (names(cur.limits) == "max") {
            sel.idx <- which(values <= cur.limits)
            val <- values[sel.idx]
            if (length(val) != 0) {
                new.setting["max"] <- unname(quantile(val, relax.speed))
            } else {
                # no values within the current limits -> set constraint to the limit
                return(cur.limits)
            }
            if (new.setting["max"] == cur.setting["max"]) {
                new.setting["max"] <- max(val, na.rm = TRUE)
            }
        } else {
            # min constraint
            sel.idx <- which(values >= cur.limits)
            val <- values[sel.idx]
            if (length(val) != 0) {
                new.setting["min"] <- unname(quantile(val, 1 - relax.speed))
            } else {
                return(cur.limits)  # nothing to relax 
            }
            if (new.setting["min"] == cur.setting["min"]) {
                new.setting["min"] <- min(val, na.rm = TRUE)
            }
        }
    }

    return(new.setting)
}
#' Filtering of Primers
#'
#' Filters a primer set during the optimization procedure.
#' 
#' @param primer.df Primer data frame.
#' @param sample Name of the current template sample.
#' @param template.df Template data frame.
#' @param settings List with settings for the constraints to be used for filtering.
#' @param mode.directionality Primer direction.
#' @param required.cvg Required ratio of covered templates.
#' If \code{required.cvg} is set to 0, the constraints are not relaxed.
#' @param target.temps Target melting temperature of the primers in Celsius.
#' This argument is only required if we try to match the melting temperatures of another primer set,
#' e.g. when first optimizing forward and then optimizing reverse primers.
#' @param results.loc Path to a directory where the results should be written.
#' @return The filtered primer data frame with respect to \code{required.cvg}.
#' @keywords internal
filter.primer.set.opti <- function(primer.df, sample, template.df, settings, 
                                   mode.directionality, required.cvg, results.loc, 
                                   target.temps) {
    
    if (length(primer.df) == 0) {
        return(NULL)
    }
    filtered.result <- cascaded.filter(primer.df, template.df, settings, mode.directionality, 
        required.cvg, target.temps, results.loc = results.loc)
    if (length(filtered.result) == 0) {
        # no feasible solution
        return(NULL)
    }
    if (length(results.loc) != 0) {
        # store data
        filtered.df <- filtered.result$data
        filtered.stats <- filtered.result$stats
        stats.relax <- filtered.result$stats_relax
        template.df <- update_template_cvg(template.df, filtered.df, mode.directionality)
        write.csv(template.df, file.path(results.loc, paste(sample, "_filtered_templates.csv", 
            sep = "")), row.names = FALSE)
        write.csv(filtered.df, file.path(results.loc, paste(sample, "_filtered_primers.csv", 
            sep = "")), row.names = FALSE)
        write.csv(filtered.stats, file.path(results.loc, paste(sample, "_filtered_primers_stats.csv", 
            sep = "")), row.names = FALSE)
        write.csv(stats.relax, file.path(results.loc, paste(sample, "_filtered_primers_stats_relax.csv", 
            sep = "")), row.names = FALSE)
        excluded.df <- filtered.result$excluded
        write.csv(excluded.df, file.path(results.loc, paste(sample, "_excluded_primers.csv", 
            sep = "")), row.names = FALSE)
    }
    return(filtered.result)
}

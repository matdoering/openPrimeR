############
# Evaluation functions for determining whether constraints are fulfilled
#############

#' Apply Constraints to a List.
#'
#' Checks whether the input values are within the specified limits.
#'
#' Applies a constraint to every element in a vector of comma separated strings.
#' Applied when filtering covered seqs according to primer efficiency.
#' 
#' @param gc.ratio.fw Forward values (comma-separated strings).
#' @param gc.ratio.rev Reverse values (comma-separated strings).
#' @param min.GC Minimal allowed value.
#' @param max.GC Maximal allowed value.
#' @param fw.idx Indices of forward values to consider.
#' @param rev.idx Indices of reverse values to consider.
#' @param mode.directionality Direction of primers
#' @return Data frame with \code{TRUE} for values fulfilling the constraints, 
#' \code{FALSE} otherwise.
#' @keywords internal
apply.constraint.list <- function(gc.ratio.fw, gc.ratio.rev, min.GC, max.GC, fw.idx, 
    rev.idx, mode.directionality = c("fw", "rev", "both")) {

    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    if (is.character(gc.ratio.fw)) {
        gc.ratio.fw <- lapply(gc.ratio.fw, function(x) as.numeric(unlist(strsplit(x, 
            split = ","))))
        gc.ratio.rev <- lapply(gc.ratio.rev, function(x) as.numeric(unlist(strsplit(x, 
            split = ","))))
    }
    if (length(min.GC) != 0 && is.na(min.GC)) {
        min.GC <- numeric(0)
    }
    if (length(max.GC) != 0 && is.na(max.GC)) {
        max.GC <- numeric(0)
    }
    if (length(min.GC) != 0 && length(max.GC) != 0) {
        s <- min.GC
        e <- max.GC
    } else if (length(min.GC) == 0) {
        s <- -Inf
        e <- max.GC
    } else {
        s <- min.GC
        e <- Inf
    }
    if (length(fw.idx) != 0) {
        b.fw <- lapply(gc.ratio.fw[fw.idx], function(x) x >= s & x <= e)
    } else {
        b.fw <- lapply(gc.ratio.fw, function(x) x >= s & x <= e)
    }
    if (length(rev.idx) != 0) {
        b.rev <- lapply(gc.ratio.rev[rev.idx], function(x) x >= s & x <= e)
    } else {
        b.rev <- lapply(gc.ratio.rev, function(x) x >= s & x <= e)
        
    }
    res.fw <- vector("list", length(gc.ratio.fw))
    res.rev <- vector("list", length(gc.ratio.rev))
    if (length(fw.idx) != 0) {
        res.fw[fw.idx] <- b.fw
    } else {
        res.fw <- b.fw
    }
    if (length(rev.idx) != 0) {
        res.rev[rev.idx] <- b.rev
    } else {
        res.rev <- b.rev
    }
    res <- data.frame(Forward = rep("", length(gc.ratio.fw)), Reverse = rep("", length(gc.ratio.fw)))
    if (mode.directionality == "both") {
        res$Forward <- string.rep(res.fw)
        res$Reverse <- string.rep(res.rev)
    } else if (mode.directionality == "fw") {
        res$Forward <- string.rep(res.fw)
    } else {
        res$Reverse <- string.rep(res.rev)
    }
    return(res)
}
#' Application of Constraints
#'
#' Checks whether the input values are within the specified limits.
#'
#' @param gc.ratio.fw Forward values.
#' @param gc.ratio.rev Reverse values.
#' @param min.GC Minimal allowed value.
#' @param max.GC Maximal allowed value.
#' @param fw.idx Indices of forward values to consider.
#' @param rev.idx Indices of reverse values to consider.
#' @param mode.directionality Direction of primers
#' @return Data frame with \code{TRUE} for values fulfilling the constraints, 
#' \code{FALSE} otherwise.  Also returns \code{FALSE} if a data point is not available.
#' @keywords internal
apply.constraint <- function(gc.ratio.fw, gc.ratio.rev, min.GC, max.GC, fw.idx, rev.idx, 
    mode.directionality = c("fw", "rev", "both")) {
    if (length(min.GC) == 0 && length(max.GC) == 0) {
        return(data.frame(rep(TRUE, length(gc.ratio.fw))))
    }
    if (length(min.GC) != 0 && is.na(min.GC)) {
        min.GC <- numeric(0)
    }
    if (length(max.GC) != 0 && is.na(max.GC)) {
        max.GC <- numeric(0)
    }
    if (length(min.GC) != 0 && length(max.GC) != 0) {
        s <- min.GC
        e <- max.GC
    } else if (length(min.GC) == 0) {
        s <- -Inf
        e <- max.GC
    } else {
        s <- min.GC
        e <- Inf
    }
    result <- rep(TRUE, length(gc.ratio.fw))
    r.fw <- gc.ratio.fw[fw.idx] >= s & gc.ratio.fw[fw.idx] <= e
    res.fw <- rep(TRUE, length(gc.ratio.fw))
    if (length(fw.idx) != 0) {
        # only take constraint values if input was given (!= '') as determined by fw.idx
        res.fw[fw.idx] <- r.fw
    }
    r.rev <- gc.ratio.rev[rev.idx] >= s & gc.ratio.rev[rev.idx] <= e
    res.rev <- rep(TRUE, length(gc.ratio.rev))
    if (length(rev.idx) != 0) {
        res.rev[rev.idx] <- r.rev
    }
    # ensure that not avaialable values are FALSE
    res.rev[is.na(res.rev)] <- FALSE
    res.fw [is.na(res.fw)] <- FALSE
    if (mode.directionality == "both") {
        res <- data.frame(Forward = res.fw, Reverse = res.rev)
    } else if (mode.directionality == "fw") {
        res <- data.frame(Forward = res.fw)
    } else {
        res <- data.frame(Reverse = res.rev)
    }
    return(res)
}
string.rep <- function(values) {
    str <- unlist(lapply(values, function(x) {
        if (length(x) != 0 && all(is.na(x))) { 
            ""
        } else {
            paste0(x, collapse = ",")
        }
        }))
    if (length(str) == 0) {
        # ensure we output a string
        str <- ""
    }
    return(str)
}
#' Computation of Constraints.
#'
#' Computes the specified constraints for the input primers.
#'
#' @param primer.df Primer data frame.
#' @param mode.directionality Primer directionality.
#' @param template.df Template data frame.
#' @param settings A \code{DesignSettings} object.
#' @param active.constraints Strings giving the constraints that are to be computed.
#' @param no.structures Whether dimer structures shall be computed.
#' @param for.shiny Whether to format output for HTML.
#' @param updateProgress Progress callback function for shiny.
#' @return A data frame with columns for each constraint in \code{active.constraints}.
#' @keywords internal
compute.constraints <- function(primer.df, mode.directionality = c("fw", "rev", "both"), template.df, settings,
    active.constraints = c("primer_coverage", "primer_length", "primer_specificity",
                            "gc_clamp", "gc_ratio", "no_runs", "no_repeats", "self_dimerization", 
                            "cross_dimerization", "melting_temp_range", 
                            "melting_temp_diff", "secondary_structure", 
                            # coverage conditions:
                            "primer_efficiency", "annealing_DeltaG",
                            "stop_codon", "terminal_mismatch_pos",
                            "substitution", "hexamer_coverage",
                            "coverage_model",
                            # off-target cvg constraints
                            "off_primer_efficiency", "off_annealing_DeltaG",
                            "off_coverage_model"),
    no.structures = FALSE, for.shiny = FALSE, updateProgress = NULL) {
    if (length(active.constraints) == 0) {
        return(NULL)
    }
    ###################
    # retrieve settings
    primer_conc <- PCR(settings)$primer_concentration
    template_conc <- PCR(settings)$template_concentration
    na_salt_conc <- PCR(settings)$Na_concentration
    mg_salt_conc <- PCR(settings)$Mg_concentration
    k_salt_conc <- PCR(settings)$K_concentration
    tris_salt_conc <- PCR(settings)$Tris_concentration
    annealing.temp <- PCR(settings)$annealing_temp
    if (length(annealing.temp) == 1) {
        # ensure that every primer has an individual annealing temperature
        annealing.temp <- rep(annealing.temp, nrow(primer.df))
    }
    taqEfficiency <- PCR(settings)$use_taq_polymerase
    allowed.mismatches <- conOptions(settings)$allowed_mismatches
    allowed.other.binding.ratio <- conOptions(settings)$allowed_other_binding_ratio
    allowed.region.definition <- conOptions(settings)$allowed_region_definition 
    ##################
    if (length(active.constraints) != 0) { # don't match if there are no constraints
        active.constraints <- match.arg(active.constraints, several.ok = TRUE)
    }
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    # Note: only place values with different names than original constraints arising from active constraints directly here (e.g.
    # with the same name as constraints. others will be removed.  store all
    # constraint results in 'constraint.values'
    constraint.values <- data.frame(#mean_primer_efficiency = numeric(nrow(primer.df)),
        primer_length_fw = numeric(nrow(primer.df)), primer_length_rev = numeric(nrow(primer.df)), 
        gc_clamp_fw = numeric(nrow(primer.df)), gc_clamp_rev = numeric(nrow(primer.df)), 
        gc_ratio_fw = numeric(nrow(primer.df)), gc_ratio_rev = numeric(nrow(primer.df)), 
        no_repeats_fw = numeric(nrow(primer.df)), no_repeats_rev = numeric(nrow(primer.df)), 
        no_runs_fw = numeric(nrow(primer.df)), no_runs_rev = numeric(nrow(primer.df)))
    original.cols <- colnames(constraint.values)
    melting.temps <- NULL  # required for specificity
    if ("melting_temp_range" %in% active.constraints || "melting_temp_diff" %in% active.constraints) {
        # melting temperature
        if (is.function(updateProgress)) {
            detail <- "melting temp"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        # melting temperature T_m:
        melting.temps <- compute.melting.temps(primer.df, primer_conc, 
                na_salt_conc, mg_salt_conc, k_salt_conc, 
                tris_salt_conc, mode.directionality)
        # update melting temp
        primer.df$melting_temp <- melting.temps$melting_temp
        # prevent an error due to cbind with data frame is empty
        if (length(constraint.values) == 0) {
            constraint.values <- melting.temps
        } else {
            constraint.values <- cbind(constraint.values, melting.temps)
        }
    }
    # determine annealing temp if required but not available
    annealing.required <- any(active.constraints %in% c("secondary_structure", "primer_efficiency", 
                                                        "annealing_DeltaG", "off_primer_efficiency",
                                                        "off_annealing_DeltaG",
                                                        "self_dimerization", "cross_dimerization"))
    if (length(annealing.temp) == 0 && annealing.required) {
        annealing.temp <- compute_annealing_temp(primer.df, mode.directionality, 
                  template.df, na_salt_conc, mg_salt_conc, k_salt_conc, 
                  tris_salt_conc, primer_conc)
    }
    cvg.required <- FALSE
    if ("primer_coverage" %in% active.constraints || "primer_specificity" %in% active.constraints && !"primer_specificity" %in% colnames(primer.df)) {
        # only trigger cvg computation for specificity purposes if coverage is explicitly required by constraints, otherwise use existing entries.
        cvg.required <- TRUE
    }
    if (cvg.required) {
        # covered templates by primers
        if (is.function(updateProgress)) {
            detail <- "primer coverage"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        # evaluate primer cvg first to have coverage data all constraints (specificity)
        cvg.df <- evaluate.primer.cvg(template.df, primer.df, mode.directionality, settings, updateProgress)
        # update primer df with coverage (might be needed for other constraints).
        primer.df <- update.constraint.values(primer.df, cvg.df)
        constraint.values <- cbind(constraint.values, cvg.df)
    }
    if ("primer_length" %in% active.constraints) {
        # length of primers
        if (is.function(updateProgress)) {
            detail <- "primer length"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        # should already be computed in the input, no computation necessary
        if (!"primer_length_fw" %in% colnames(primer.df)) {
            stop("WARNING: primer len not computed already")
        }
        constraint.values$primer_length_fw <- primer.df$primer_length_fw
        constraint.values$primer_length_rev <- primer.df$primer_length_rev
    }
    if ("gc_clamp" %in% active.constraints) {
        # count of gc at 3' end
        if (is.function(updateProgress)) {
            detail <- "gc clamp"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        gc.count.fw <- evaluate.GC.clamp(primer.df$Forward)
        gc.count.rev <- evaluate.GC.clamp(primer.df$Reverse)
        constraint.values$gc_clamp_fw <- gc.count.fw
        constraint.values$gc_clamp_rev <- gc.count.rev
    }
    if ("gc_ratio" %in% active.constraints) {
        # ratio of gc
        if (is.function(updateProgress)) {
            detail <- "gc ratio"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        gc.ratio.fw <- compute.gc.ratio(primer.df$Forward)
        gc.ratio.rev <- compute.gc.ratio(primer.df$Reverse)
        constraint.values$gc_ratio_fw <- gc.ratio.fw
        constraint.values$gc_ratio_rev <- gc.ratio.rev
    }
    if ("no_repeats" %in% active.constraints) {
        # number of repeats
        if (is.function(updateProgress)) {
            detail <- "no repeats"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        no.repeats.fw <- nbr.of.repeats(primer.df$Forward)
        no.repeats.rev <- nbr.of.repeats(primer.df$Reverse)
        constraint.values$no_repeats_fw <- no.repeats.fw
        constraint.values$no_repeats_rev
    }
    if ("no_runs" %in% active.constraints) {
        # number of runs
        if (is.function(updateProgress)) {
            detail <- "no runs"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        nbr.runs.fw <- nbr.of.runs(primer.df$Forward)
        nbr.runs.rev <- nbr.of.runs(primer.df$Reverse)
        constraint.values$no_runs_fw <- nbr.runs.fw
        constraint.values$no_runs_rev <- nbr.runs.rev
    }
    #######
    # COVERAGE CONSTRAINTS: terminal mismatches, stop codons, efficiency, annealingDeltaG
    #######
    if (any(c("stop_codon", "substitution") %in% active.constraints)) {
        if (is.function(updateProgress)) {
            detail <- "Identifying mutations"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        #######################
        # Stop codons
        ###################
        mutation.types <- c("stop_codon", "substitution")[c("stop_codon", "substitution") %in% active.constraints]
        mutation.check <- mismatch.mutation.check(primer.df, template.df, mutation.types)
        out.strings <- lapply(mutation.check, function(x) 
            if (length(x) == 0) { 
                ""
            } else {
                unlist(lapply(seq_len(ncol(x)), function(y) paste(x[, y], collapse =",")))
        })
        if ("stop_codon" %in% mutation.types) {
            stop.str <- sapply(out.strings, function(x) x[1])
            constraint.values$stop_codon <- stop.str
        }
        if ("substitution" %in% mutation.types) {
            sub.str <- sapply(out.strings, function(x) x[2]) 
            constraint.values$substitution <- sub.str
        }
    }
    if ("terminal_mismatch_pos" %in% active.constraints) {
        ###########################
        # Identify 3' mismatches
        ############################
        if (is.function(updateProgress)) {
            detail <- "Terminal mismatch identification"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        mm.check <- check.3prime.mismatches(template.df, primer.df, mode.directionality)
        mm.str <- string.rep(mm.check)
        constraint.values$terminal_mismatch_pos <- mm.str
    }
    if ("hexamer_coverage" %in% active.constraints && length(conOptions(settings)$hexamer_coverage) != 0 && conOptions(settings)$hexamer_coverage == "active") {
        # only computed when active because presence of hexamer_coverage is used for ignoring the filtering later! 
        ###########################
        # Identify 3' hexamers with 100% complementarity to template
        ############################
        if (is.function(updateProgress)) {
            detail <- "3' Hexamer detection"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        hex.check <- check.3prime.hexamers(template.df, primer.df, mode.directionality)
        hex.str <- string.rep(hex.check)
        constraint.values$hexamer_coverage <- hex.str
    }
    if ("coverage_model" %in% active.constraints) {
        #######
        # Logistic regression model for predicting coverage
        ########
        if (is.function(updateProgress)) {
            detail <- "Coverage model"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        preds <- predict_coverage(primer.df, template.df, settings,
                                  updateProgress = updateProgress)
        pred.string <- string.rep(preds)
        constraint.values$coverage_model <- pred.string
    }
    if ("off_coverage_model" %in% active.constraints) {
        #######
        # Logistic regression model for predicting coverage
        ########
        if (is.function(updateProgress)) {
            detail <- "Coverage model"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        preds <- predict_coverage(primer.df, template.df, settings, mode = "off_target",
                                  updateProgress = updateProgress)
        pred.string <- string.rep(preds)
        constraint.values$off_coverage_model <- pred.string
    }
    if ("primer_efficiency" %in% active.constraints) {
        ############
        # Thermodynamic model of annealing + 3' terminal mismatch model (for taq) 
        #############
        if (is.function(updateProgress)) {
            detail <- "Primer efficiency"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        eff_taq <- compute.primer.efficiencies(primer.df, template.df, 
                            annealing.temp, taqEfficiency = TRUE, 
                            PCR(settings)$primer_concentration,
                            PCR(settings)$Na_concentration,
                            PCR(settings)$Mg_concentration,
                            PCR(settings)$K_concentration,
                            PCR(settings)$Tris_concentration,
                            eff.only = FALSE, mode = "on_target")
        eff.string <- string.rep(eff_taq)
        eff.mean <- sapply(eff_taq, function(x) ifelse(length(x) == 0, 0, mean(x)))
        constraint.values$primer_efficiency <- eff.string
        constraint.values$mean_primer_efficiency <- eff.mean
    }
    if ("off_primer_efficiency" %in% active.constraints) {
        ############
        # Off-target binding events: Thermodynamic model of annealing + 3' terminal mismatch model (for taq) 
        #############
        if (is.function(updateProgress)) {
            detail <- "Primer efficiency"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        eff_off <- compute.primer.efficiencies(primer.df, template.df, 
                            annealing.temp, taqEfficiency = TRUE, 
                            PCR(settings)$primer_concentration,
                            PCR(settings)$Na_concentration,
                            PCR(settings)$Mg_concentration,
                            PCR(settings)$K_concentration,
                            PCR(settings)$Tris_concentration,
                            eff.only = FALSE, mode = "off_target")
        eff.string <- string.rep(eff_off)
        constraint.values$off_primer_efficiency <- eff.string
    }
    if ("annealing_DeltaG" %in% active.constraints) {
        ######################
        # Duplex binding energies
        ########################
        if (is.function(updateProgress)) {
            detail <- "Free energy of annealing"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        deltaG_ABs <- get.duplex.energies(primer.df, template.df, annealing.temp, settings, mode = "on_target")
        deltaG.string <- string.rep(deltaG_ABs)
        deltaG.mean <- unlist(lapply(deltaG_ABs, function(x) ifelse(length(x) == 0, 0, mean(x))))
        constraint.values$annealing_DeltaG <- deltaG.string
        constraint.values$mean_annealing_DeltaG <- deltaG.mean
    }
     if ("off_annealing_DeltaG" %in% active.constraints) {
        ######################
        # Off-target duplex binding energies
        ########################
        if (is.function(updateProgress)) {
            detail <- "Off-targets: free energy of annealing"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        deltaG_ABs <- get.duplex.energies(primer.df, template.df, annealing.temp, settings, mode = "off_target") # TODO: off target nbr is wrong ..
        deltaG.string <- string.rep(deltaG_ABs)
        constraint.values$off_annealing_DeltaG <- deltaG.string
    }
    # get secondary structures
    if ("secondary_structure" %in% active.constraints) {
        # primer secondary structures
        if (is.function(updateProgress)) {
            detail <- "secondary structure"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        structs <- compute.secondary.structures(primer.df, mode.directionality, annealing.temp)
        if (length(constraint.values) == 0) {
            constraint.values <- structs
        } else {
            constraint.values <- cbind(constraint.values, structs)
        }
    }
    if ("primer_specificity" %in% active.constraints) {
        # primer specificity
        if (length(primer.df$primer_specificity) != 0) {
            # look at existing data
            constraint.values$primer_specificity <- primer.df$primer_specificity
        } 
    }
    if ("self_dimerization" %in% active.constraints) {
        if (is.function(updateProgress)) {
            detail <- "self dimerization"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        self.dimer.df <- compute.all.self.dimers.frontend(primer.df, primer_conc, 
            na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, annealing.temp, for.shiny = for.shiny, no.structures = no.structures)
        if (length(self.dimer.df) != 0) {
            # we don't have any column entries when there's nothing to report ..
            if (length(constraint.values) == 0 || nrow(constraint.values) == 0) {
                constraint.values <- self.dimer.df
            } else {
                constraint.values <- cbind(constraint.values, self.dimer.df)
            }
        }
    }
    if ("cross_dimerization" %in% active.constraints) {
        if (is.function(updateProgress)) {
            detail <- "cross dimerization"
            updateProgress(1/length(active.constraints), detail, "inc")
        }
        cross.dimer.df <- compute.all.cross.dimers.frontend(primer.df, primer_conc, 
            na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, annealing.temp, for.shiny = for.shiny, no.structures = no.structures)
        if (length(cross.dimer.df) != 0) {
            # no entries in constraint.values when there's nothing to report ...
            if (length(constraint.values) == 0 || nrow(constraint.values) == 0) {
                constraint.values <- cross.dimer.df
            } else {
                constraint.values <- cbind(constraint.values, cross.dimer.df)
            }
        }
    }
    # remove all 'constraint.values' that that were not actually computed
    keep.idx.a <- which(!(colnames(constraint.values) %in% original.cols))
    keep.idx.b <- unlist(lapply(active.constraints, function(x) grep(x, colnames(constraint.values))))
    keep.idx.c <- unlist(lapply(names(cvg_constraints(settings)), function(x) grep(x, colnames(constraint.values)))) # retain cvg constraint values
    #print("to keep:")
    #print(colnames(constraint.values[keep.idx.a]))
    keep.idx <- union(keep.idx.a, c(keep.idx.b, keep.idx.c))
    #print(colnames(constraint.values)[keep.idx])
    constraint.values <- constraint.values[, keep.idx, drop = FALSE]
    return(constraint.values)
}
#' Check Constraints on Templates
#'
#' Transforms the comma-separated input strings to a boolean representation.
#'
#' @param template.constraints Strings with comma-separated values to be turned to logical.
#' @return List with boolean values
#' @keywords internal
check.template.constraints <- function(template.constraints) {
    # evaluate list constraints (e.g. primer efficiency)
    if (length(template.constraints) == 0) {
        return(NULL)
    }
    all.res <- vector("list", length(template.constraints))
    # message('check.template.constraints:') message(template.constraints) for each
    # template constraint
    for (i in seq_along(template.constraints)) {
        data <- template.constraints[[i]]
        res <- lapply(seq_len(nrow(data)), function(y) {
            # each x is two lists, one for fw, one for rev
            x <- data[y, ]
            # message(x)
            s1 <- as.logical(strsplit(as.character(x[[1]]), split = ",")[[1]])
            s2 <- as.logical(strsplit(as.character(x[[2]]), split = ",")[[1]])
            if (length(s1) == 0) {
                res <- s2
            } else if (length(s2) == 0) {
                res <- s1
            } else {
                res <- unlist(lapply(seq_along(s1), function(y) s1[y] && s2[y]))
            }
        })
        all.res[[i]] <- res
    }
    names(all.res) <- names(template.constraints)
    return(all.res)
}
#' Merge Template Decisions.
#' 
#' Merges the results for multiple template evaluations.
#'
#' @param eval.t List with evaluated template constraints
#' @return List with merged boolean decisions.
#' @keywords internal
merge.template.decisions <- function(eval.t) {
    # get merged result for all template decisions
    final.res <- NULL
    if (length(eval.t) != 0) {
        final.res <- vector("list", length(eval.t[[1]]))
    }
    for (i in seq_along(final.res)) {
        # for every primer
        template.decisions <- rep(TRUE, length(eval.t[[1]][[i]]))
        for (j in seq_along(eval.t)) {
            # for every constraint
            cur.decisions <- eval.t[[j]][[i]]
            f.idx <- which(!cur.decisions)
            if (length(f.idx) != 0) {
                template.decisions[f.idx] <- FALSE
            }
        }
        final.res[[i]] <- template.decisions
    }
    return(final.res)
}
#' Check for Relaxation
#'
#' Determines whether constraints where relaxed or not.
#'
#' @param used.constraints The constraints that were used during the optimization.
#' @param input.constraints The user-specified constraints.
#' @return If \code{input.constraints} was relaxed TRUE is returned, otherwise FALSE.
#' @keywords internal
were.constraints.relaxed <- function(used.constraints, input.constraints) {
    # determines whether used constraints differ from input constraints
    
    m <- match(names(used.constraints), names(input.constraints))
    was.relaxed <- any(sapply(seq_along(m), function(x) any(used.constraints[[x]] != 
        input.constraints[[m[x]]])))
    return(was.relaxed)
}
#' Evaluation of Template Constraints.
#'
#' Evaluates the input template constraints.
#'
#' @param constraint.values Data frame with template constraints
#' @param constraint.settings List specifying the allowed values for constraint evaluation.
#' @param active.constraints Strings specifying the constraints to check.
#' @param mode.directionality Direction of primers.
#' @return List indicating which template constraints where fulfilled or not (TRUE/FALSE).
#' @keywords internal
evaluate.template.constraints <- function(constraint.values, constraint.settings, 
                                          active.constraints, mode.directionality = c("fw", "rev", "both")) {

    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    template.constraints <- list()  # constraints on template-specific level (lists)
    if ("coverage_model" %in% active.constraints) {
        cvg.model <- apply.constraint.list(constraint.values$coverage_model, 
            constraint.values$coverage_model, constraint.settings[["coverage_model"]]["min"], 
            constraint.settings[["coverage_model"]]["max"], NULL, NULL, mode.directionality)
        template.constraints[["T_EVAL_coverage_model"]] <- cvg.model
    }
    if ("off_coverage_model" %in% active.constraints) {
        model.FPR <- apply.constraint.list(constraint.values$off_coverage_model, 
            constraint.values$off_coverage_model, constraint.settings[["off_coverage_model"]]["min"],
            constraint.settings[["off_coverage_model"]]["max"], NULL, NULL, mode.directionality)
        template.constraints[["T_off_EVAL_coverage_model"]] <- model.FPR
    }

    if ("primer_efficiency" %in% active.constraints) {
        primer.eff <- apply.constraint.list(constraint.values$primer_efficiency, 
            constraint.values$primer_efficiency, constraint.settings[["primer_efficiency"]]["min"], 
            constraint.settings[["primer_efficiency"]]["max"], NULL, NULL, mode.directionality)
        template.constraints[["T_EVAL_primer_efficiency"]] <- primer.eff
    }
    if ("off_primer_efficiency" %in% active.constraints) { # off-target events
        primer.eff <- apply.constraint.list(constraint.values$off_primer_efficiency, 
            constraint.values$off_primer_efficiency, constraint.settings[["off_primer_efficiency"]]["min"], 
            constraint.settings[["off_primer_efficiency"]]["max"], NULL, NULL, mode.directionality)
        template.constraints[["T_off_EVAL_primer_efficiency"]] <- primer.eff
    }
    if ("annealing_DeltaG" %in% active.constraints) {
        annealing.deltaG <- apply.constraint.list(constraint.values$annealing_DeltaG, 
            constraint.values$annealing_DeltaG, constraint.settings[["annealing_DeltaG"]]["min"],
            constraint.settings[["annealing_DeltaG"]]["max"], NULL, NULL, mode.directionality)
        template.constraints[["T_EVAL_annealing_DeltaG"]] <- annealing.deltaG
    }
    if ("off_annealing_DeltaG" %in% active.constraints) {
        annealing.deltaG <- apply.constraint.list(constraint.values$off_annealing_DeltaG, 
            constraint.values$off_annealing_DeltaG, constraint.settings[["off_annealing_DeltaG"]]["min"],
            constraint.settings[["off_annealing_DeltaG"]]["max"], NULL, NULL, mode.directionality)
        template.constraints[["T_off_EVAL_annealing_DeltaG"]] <- annealing.deltaG
    }
    if ("stop_codon" %in% active.constraints) {
        stop.sel <- apply.constraint.list(constraint.values$stop_codon, constraint.values$stop_codon, 
                                        constraint.settings[["stop_codon"]]["min"], constraint.settings[["stop_codon"]]["max"],
                                        NULL, NULL, mode.directionality)
        template.constraints[["T_EVAL_stop_codon"]] <- stop.sel
    }
    if ("substitution" %in% active.constraints) {
        sub.sel <- apply.constraint.list(constraint.values$substitution, constraint.values$substitution, 
                                        constraint.settings[["substitution"]]["min"], constraint.settings[["substitution"]]["max"],
                                        NULL, NULL, mode.directionality)
        template.constraints[["T_EVAL_substitution"]] <- sub.sel
    }

    if ("terminal_mismatch_pos" %in% active.constraints) {
        # 3' Mismatch filtering
        mm.sel <- apply.constraint.list(constraint.values$terminal_mismatch_pos, constraint.values$terminal_mismatch_pos, 
                                        constraint.settings[["terminal_mismatch_pos"]]["min"], constraint.settings[["terminal_mismatch_pos"]]["max"],
                                        NULL, NULL, mode.directionality)

        template.constraints[["T_EVAL_terminal_mismatch_pos"]] <- mm.sel
    } 
    template.con.result <- check.template.constraints(template.constraints)
    return(template.con.result)
}

#' Update of Primer Constraints.
#'
#' Updates the input primer data frame with the computed constraint values.
#'
#' @param constraint.df Primer data frame.
#' @param constraint.values Data frame with computed constraint values.
#' @return A primer data frame with updated columns.
#' @keywords internal
update.constraint.values <- function(constraint.df, constraint.values) {
    # update/insert computed constraint.values into constraint.df
    if (length(constraint.values) != 0) {
        idx.replace <- which(colnames(constraint.values) %in% colnames(constraint.df))
        idx.add <- which(!(colnames(constraint.values) %in% colnames(constraint.df)))
        if (length(idx.replace) != 0) {
            constraint.df[, colnames(constraint.values)[idx.replace]] <- constraint.values[, 
                idx.replace]
        }
        if (length(idx.add) != 0) {
            constraint.df[, colnames(constraint.values)[idx.add]] <- constraint.values[, 
                idx.add]
        }
    }
    return(constraint.df)
}
#' Evaluation of Primer Constraints.
#'
#' Determines whether a set of primers
#' fulfills the constraints on the properties of the primers.
#' 
#' When the optional argument
#' \code{active.constraints} is supplied, only a subset of the constraints
#' provided in \code{settings} is evaluated. Only constraints that
#' are defined in \code{settings} can be computed. For a detailed
#' description of all possible constraints and their options, please
#' consider the \code{\link{ConstraintSettings}} documentation.
#'
#' @param primer.df A \code{Primers} object containing the primers
#' whose properties are to be checked.
#' @param template.df A \code{Templates} object containing the 
#' template sequences corresponding to \code{primer.df}.
#' @param settings A \code{DesignSettings} object containing the 
#' constraints that are to be evaluated.
#' @param active.constraints A subset of the constraint identifiers 
#' provided by \code{settings} that are to be checked
#' for fulfillment. By default \code{active.constraints} is \code{NULL} such that
#' all constraints found in \code{settings} are evaluated. Otherwise,
#' only the constraints specified via \code{active.constraints} 
#' that are available in \code{settings} are considered.
#' @param to.compute.constraints Constraints that are to be computed.
#' By default, \code{to.compute.constraints} is set to \code{NULL} such that
#' all \code{active.constraints} are computed. If \code{to.compute.constraints}
#' is a subset of \code{active.constraints}, all constraints specified
#' via \code{active.constraints} are evaluated for fulfillment,
#' but only the constraints in \code{to.compute.constraints} are newly calculated.
#' @param for.shiny Whether the output of the function shall be
#' formatted as HTML. The default setting is \code{FALSE}.
#' @param updateProgress Progress callback function for shiny. The defaut is
#' \code{NULL} meaning that no progress is monitored via the Shiny interface.
#' @return A \code{Primers} object that is augmented
#' with columns indicating the results for each evaluated constraint.
#' The \code{constraints_passed} column indicates whether all \code{active.constraints} were fulfilled.
#' The \code{EVAL_*} columns indicate the fulfillment of primer-specific constraints.
#' The \code{T_EVAL_*} columns indicate the fulfillment of template-specific
#' (e.g. coverage-based) constraints.
#' For the coverage computations, columns prefixed by \code{Basic_},
#' indicate the results from string matching, while all other results
#' (e.g. \code{primer_coverage}) indicate the expected coverage
#' after applying the coverage constraints specified in \code{settings}.
#' Columns prefixed by \code{Off_} indicate off-target binding results.
#' @note Please note that some constraints can only be computed if additional software is installed,
#' please see the documentation of
#' \code{\link{DesignSettings}} for an overview.
#' @export
#' @keywords Primers
#' @family primer functions
#' @examples
#' data(Ippolito)
#' settings.xml <- system.file("extdata", "settings", 
#'                  "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
#' settings <- read_settings(settings.xml)
#' # Check all constraints found in 'settings':
#' constraint.df <- check_constraints(primer.df, template.df, 
#'                      settings, active.constraints = names(constraints(settings)))
#' # Summarize the evaluation results
#' summary(constraint.df)
check_constraints <- function(primer.df, template.df, settings,
    active.constraints = names(constraints(settings)),
    to.compute.constraints = active.constraints, for.shiny = FALSE, updateProgress = NULL) {
    if (length(settings) == 0 || !is(settings, "DesignSettings")) {
        stop("check_constraints: Please provide a 'DesignSettings' object for 'settings'.")
    }
    if (!is(primer.df, "Primers")) {
        stop("Please input a 'Primers' object for 'primer.df'.")
    }
    if (!is(template.df, "Templates")) {
        stop("Please input a valid 'Templates' object for 'template.df'.")
    }
    mode.directionality <- get.analysis.mode(primer.df)
    # only allow active.constraints for which a setting exists
    ok.constraints <- active.constraints %in% names(constraints(settings))
    bad.constraints <- which(!ok.constraints)
    if (length(bad.constraints) != 0) {
        warning("The following constraints were not found in 'settings' and excluded: \n",
                paste0(active.constraints[bad.constraints], collapse = ", ")) 
    }
    active.constraints <- active.constraints[ok.constraints]
    if (is.null(to.compute.constraints)) {
        # if NULL, compute all constraints
        to.compute.constraints <- active.constraints
    }
    to.compute.constraints <- to.compute.constraints[to.compute.constraints %in% active.constraints]
    if (length(active.constraints) == 0) {
        # nothing to check/compute
        warning("No 'active.constraints' available: no constraints checked!")
        primer.df$constraints_passed <- rep(TRUE, nrow(primer.df))
        return(primer.df)
    }
    PCR.settings <- PCR(settings)
    other.settings <- conOptions(settings)
    constraint.values <- compute.constraints(primer.df, mode.directionality, template.df, settings,
        to.compute.constraints, for.shiny = for.shiny, 
        updateProgress = updateProgress)
    constraint.df <- update.constraint.values(primer.df, constraint.values)  # update constraint.df with computed values
    eval <- eval.constraints(constraint.df, constraints(settings), active.constraints, 
        mode.directionality, primer.df)
    idx.replace <- which(colnames(eval) %in% colnames(constraint.df))
    idx.add <- which(!(colnames(eval) %in% colnames(constraint.df)))
    if (length(idx.replace) != 0) {
        constraint.df[, colnames(eval)[idx.replace]] <- eval[, idx.replace]
    }
    if (length(idx.add) != 0) {
        constraint.df[, colnames(eval)[idx.add]] <- eval[, idx.add]
    }
    if (length(eval) != 0) {
        not.fulfilled.idx <- which(apply(eval, 1, function(x) any(!x, na.rm = TRUE)))
    } else {
        not.fulfilled.idx <- NULL  # everything fulfilled if there's nothing to check ..
    }
    constraint.pass <- rep(TRUE, nrow(constraint.df))
    constraint.pass[not.fulfilled.idx] <- FALSE
    constraint.df$constraints_passed <- constraint.pass
    return(constraint.df)
}
#' Evaluation of Coverage Constraints.
#'
#' Computes the biochemical properties specified in the
#' \code{settings} object and determines whether the primers
#' fulfill the required constraints. 
#' 
#' @param primer.df A \code{Primers} object containing the primers
#' to be checked.
#' @param template.df A \code{Templates} object containing the 
#' template sequences corresponding to the primers.
#' @param settings A \code{DesignSettings} object containing the 
#' coverage constraints to be checked and their settings.
#' @param active.constraints Identifiers of constraints that are to be checked.
#' @param to.compute.constraints Constraints that are to be computed.
#' @param for.shiny Whether to format output for HTML.
#' @param updateProgress Progress callback function for shiny.
#' @return A \code{Primers} object with with columns for each constraint in \code{active.constraints}. 
#' @note Please note that some constraints can only be computed if additional software is installed,
#' please see \code{\link{DesignSettings-class}} for an overview.
#' @keywords internal
check_cvg_constraints <- function(primer.df, template.df, settings,
                                  active.constraints = names(cvg_constraints(settings)),
    to.compute.constraints = active.constraints, for.shiny = FALSE, updateProgress = NULL) {
    
    if (length(settings) == 0 || !is(settings, "DesignSettings")) {
        stop("check_constraints: Please check the input settings")
    }
    if (!is(primer.df, "Primers")) {
        stop("Please input a primer data frame.")
    }
    if (!is(template.df, "Templates")) {
        stop("Please input a valid template data frame.")
    }
    if (!"primer_coverage" %in% colnames(primer.df)) {
        stop("Please compute the primer coverage first before checking the coverage constraints.")
    }
    mode.directionality <- get.analysis.mode(primer.df)
    # only allow active constraints for which a setting exists
    active.constraints <- active.constraints[active.constraints %in% names(cvg_constraints(settings))]

    to.compute.constraints <- to.compute.constraints[to.compute.constraints %in% active.constraints]
    if (length(active.constraints) == 0) {
        # nothing to check/compute
        return(primer.df)
    }
    #########
    # consider on-target binding events to determine coverage
    ########
    constraint.values <- compute.constraints(primer.df, mode.directionality, template.df, settings,
        to.compute.constraints, for.shiny = for.shiny, updateProgress = updateProgress)
    constraint.df <- update.constraint.values(primer.df, constraint.values)  # update constraint.df with computed values
    eval.t <- evaluate.template.constraints(constraint.values, cvg_constraints(settings), 
                                            active.constraints, mode.directionality)
    eval.d <- lapply(eval.t, function(x) unlist(lapply(x, function(y) paste(y, collapse = ","))))
    eval.d <- do.call(cbind, eval.d)
    final.res <- merge.template.decisions(eval.t)
    idx.replace <- which(colnames(eval.d) %in% colnames(constraint.df))
    idx.add <- which(!(colnames(eval.d) %in% colnames(constraint.df)))
    if (length(idx.replace) != 0) {
        # don't use drop = TRUE here, we need individual column if we just have a single element
        constraint.df[, colnames(eval.d)[idx.replace]] <- eval.d[, idx.replace]
    }
    if (length(idx.add) != 0) {
        constraint.df[, colnames(eval.d)[idx.add]] <- eval.d[, idx.add]
    }
    # add template constraint result:
    if (length(final.res) != 0) {
        constraint.df$constraints_passed_T <- unlist(lapply(final.res, function(x) paste(x, 
            collapse = ",")))
    } 
    ##########
    # consider off-target binding constraints to update primer specificities
    ########
    # only check efficiency & deltaG constraints for off coverage events
    off.constraints <- active.constraints[active.constraints %in% c("primer_efficiency", "annealing_DeltaG", "coverage_model")]
    off.constraint.settings <- cvg_constraints(settings)[off.constraints]
    if (length(off.constraint.settings) != 0) {
        off.constraints <- paste0("off_", off.constraints) # compute constraints for off-target events 
        names(off.constraint.settings) <- off.constraints
    }
    off.constraint.vals <- compute.constraints(constraint.df, mode.directionality, template.df, settings,
                                                   off.constraints, mode.directionality)
    eval.t <- evaluate.template.constraints(off.constraint.vals, off.constraint.settings,
                                                off.constraints, mode.directionality)
    eval.d <- lapply(eval.t, function(x) unlist(lapply(x, function(y) paste(y, collapse = ","))))
    eval.d <- do.call(cbind, eval.d)
    final.res <- merge.template.decisions(eval.t)
    # add off-target info to result
    idx.replace <- which(colnames(eval.d) %in% colnames(constraint.df))
    idx.add <- which(!(colnames(eval.d) %in% colnames(constraint.df)))
    if (length(idx.replace) != 0) {
        # don't use drop = TRUE here, we need individual column if we just have a single element
        constraint.df[, colnames(eval.d)[idx.replace]] <- eval.d[, idx.replace]
    }
    if (length(idx.add) != 0) {
        constraint.df[, colnames(eval.d)[idx.add]] <- eval.d[, idx.add]
    }
    if (length(final.res) != 0) {
        # update off-target binding event values:
        constraint.df$constraints_passed_off_T <- unlist(lapply(final.res, function(x) paste(x, collapse = ",")))
    }
    return(constraint.df)
}


#' Retrive Constraint Indices.
#'
#' Gets the index of the required constraint columns in the primer data frame.
#' 
#' @param active.constraints The names of the constraints for which to find the indices in \code{constraint.df}.
#' @param constraint.df The primer data frame where the \code{active.constraints} should be found.
#' @return Indices of \code{active.constraints} in \code{constraint.df}.
#' @keywords internal
get.constraint.value.idx <- function(active.constraints, constraint.df) {
    # important: mapping from constraints to entries in constraint.df
    CONSTRAINTS.TO.VALUE.MAPPING <- list(primer_coverage = "primer_coverage", 
        primer_length = c("primer_length_fw", "primer_length_rev"), 
        primer_specificity = "primer_specificity", 
        gc_clamp = c("gc_clamp_fw", "gc_clamp_rev"), 
        gc_ratio = c("gc_ratio_fw", "gc_ratio_rev"), 
        no_runs = c("no_runs_fw", "no_runs_rev"), 
        no_repeats = c("no_repeats_fw", "no_repeats_rev"), 
        self_dimerization = "Self_Dimer_DeltaG", 
        cross_dimerization = "Cross_Dimer_DeltaG", 
        melting_temp_range = "melting_temp", 
        melting_temp_diff = "melting_temp_diff",
        secondary_structure = "Structure_deltaG", 
        # coverage constraints:
        primer_efficiency = "primer_efficiency",
        annealing_DeltaG = "annealing_DeltaG",
        stop_codon = "stop_codon",
        substitution = "substitution",
        terminal_mismatch_pos = "terminal_mismatch_pos",
        coverage_model = "coverage_model")
    if (!all(active.constraints %in% names(CONSTRAINTS.TO.VALUE.MAPPING))) {
        msg <- paste(active.constraints[!active.constraints %in% names(CONSTRAINTS.TO.VALUE.MAPPING)], collapse = ",")
        stop(paste("active.constraints: incorrect values. The following ",
                "constraints were not available: ", msg, sep = ""))
    }
    check.active <- sapply(active.constraints, function(x) any(!CONSTRAINTS.TO.VALUE.MAPPING[[x]] %in% 
        colnames(constraint.df)))
    if (any(check.active)) {
        idx <- which(check.active)
        warning(paste("Some active constraints were not contained in constraint.df according to the constraints to value mapping: ", 
            paste(active.constraints[check.active], collapse = ","), sep = ""))
    }
    value.idx <- lapply(CONSTRAINTS.TO.VALUE.MAPPING[active.constraints[!check.active]], function(x) match(x, 
        colnames(constraint.df)))
    return(value.idx)
}
#' Evaluation of Primers for Comparison
#'
#' Evaluate multiple primer sets according to the input constraint settings.
#'
#' @param primer.data List with primer data frames.
#' @param constraint.settings List with constraint.settings.
#' @return List with evaluated primer data frames.
#' @keywords internal
eval.comparison.primers <- function(primer.data, constraint.settings) {
    
    active.constraints <- names(constraint.settings)
    new.data <- primer.data
    for (i in seq_along(primer.data)) {
        constraint.df <- primer.data[[i]]
        mode.directionality <- get.analysis.mode(constraint.df)
        eval.df <- eval.constraints(constraint.df, constraint.settings, active.constraints, 
            mode.directionality, constraint.df)
        primer.df <- update.constraint.values(constraint.df, eval.df)
        new.data[[i]] <- primer.df
    }
    return(new.data)
}
#' Evaluation of Constraints'
#' 
#' Evaluates whether the given primer data frame fulfills the required conditions.
#'
#' Constraint values should be contained in \code{constraint.df}. 
#' For each constraint in \code{active.constraints}, a boolean column with the name
#' \code{EVAL_<constraint_name>} is generated, which indicates whether 
#' a primer in a given rows fulfills a constraint or not.
#'
#' @param constraint.df Primer data frame with computed constraints.
#' @param constraint.settings List with allowed values pers constraint.
#' @param active.constraints Names of constraints to be evaluated.
#' @param mode.directionality Directionality of primers
#' @param primer.df Primer data frame corresponding to \code{constraint.df}.
#' @return Augments the \code{constraint.df} data frame with evaluation columns.
#' @keywords internal
eval.constraints <- function(constraint.df, constraint.settings, active.constraints, 
    mode.directionality = c("fw", "rev", "both"), primer.df) {
    # evaluates constraints, i.e. determines whether constraints are fulfilled or
    # not. Requires that constraints have been precomputed.  Returns a data frame
    # with EVAL_<constraint_name> column, indicating whether the constraints of
    # individual primers (rows) are fulfilled or not.
   
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    # exclude constraints that shouldn't be checked
    template.constraints <- c("primer_efficiency", "annealing_DeltaG", "stop_codon", "substitution", "terminal_mismatch_pos", "hexamer_coverage", "coverage_model")
    active.constraints <- active.constraints[!active.constraints %in% template.constraints]
    if (length(active.constraints) == 0) {
        return(NULL)
    }
    # use the global variable constraints.to.value.mapping to identify the value
    # columns corresponding to the constraints:
    value.idx <- get.constraint.value.idx(active.constraints, constraint.df)  # idx of constraint values in constraint.df
    # determine if we should constrain both fw and rev, or only one direction
    fw.idx <- which(primer.df$Forward != "")
    rev.idx <- which(primer.df$Reverse != "")
    all.constraints <- vector("list", length(active.constraints))
    names(all.constraints) <- paste("EVAL_", active.constraints, sep = "")
    for (i in seq_along(value.idx)) {
        active.constraint <- names(value.idx)[i]
        #message(paste("evaluating constraint: ", active.constraint, sep = ""))
        constraint.idx <- value.idx[[i]]
        constraints <- constraint.settings[[active.constraint]]
        if (is.null(names(constraints)) || !names(constraints) %in% c("min", "max")) {
            msg <- paste("Constraint did not have a proper annotation of min/max setting: ", active.constraint, sep = "")
            stop(msg)
        }
        min.condition <- unname(constraints["min"])
        max.condition <- unname(constraints["max"])
        if (is.na(min.condition)) {
            min.condition <- NULL
        } else if (is.na(max.condition)) {
            max.condition <- NULL
        }
        if (length(constraint.idx) == 2) {
            # fw and rev values are available for current constraint
            eval <- apply.constraint(constraint.df[, constraint.idx[1]], constraint.df[, 
                constraint.idx[2]], min.condition, max.condition, fw.idx, rev.idx, 
                mode.directionality)
        } else {
            eval <- apply.constraint(constraint.df[, constraint.idx[1]], constraint.df[, 
                constraint.idx[1]], min.condition, max.condition, fw.idx, rev.idx, 
                mode.directionality)
        }
        all.constraints[[i]] <- eval
    }
    if (mode.directionality == "both") {
        all.constraints <- lapply(all.constraints, function(y) apply(y, 1, function(x) all(x, 
            na.rm = TRUE)))
    } else if (mode.directionality == "fw") {
        # message(all.constraints)
        all.constraints <- lapply(all.constraints, function(y) y[, "Forward"])
    } else {
        all.constraints <- lapply(all.constraints, function(y) y[, "Reverse"])
    }
    constraint.df <- do.call(cbind, all.constraints)
    return(constraint.df)
}


###############
# Greedy optimization of primer sets
#################

#' Reorder Primers 
#'
#' Reorders a primer set according to the IDs of primers.
#'
#' @param filtered.primers Primer data frame.
#' @param primer.ID.order new ordering of IDs in the data frame.
#' @return Reordered primer data frame.
#' @keywords internal
reorder.primer.table <- function(filtered.primers, primer.ID.order) {
    # reorders the input primers table according to the given ID order
    # filtered.primers: input primer df primer.ID.order: new order of ID column in
    # primer.df
    if (length(primer.ID.order) == 0) {
        return(filtered.primers)
    }
    filtered.primers$ID <- factor(filtered.primers$ID, levels = primer.ID.order)
    o <- order(filtered.primers$ID)
    filtered.primers <- filtered.primers[o, ]
    return(filtered.primers)
}
#' Greedy Choice
#'
#' Selects the currently best primer for Greedy primer selection.
#'
#' @param result Data frame of current optimized primer data set that is to be augmented.
#' @param primers Data frame of candidate primers for addition to \code{result}.
#' @param deltaG.cutoff Free energy cutoff for cross-dimerization.
#' @param target.temp Target annealing temperature in Celsius.
#' @param primer_conc Primer concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @return The index of a suitable primer according to Greedy selection.
#' @keywords internal
select.best.primer.idx <- function(result, primers, deltaG.cutoff, target.temp,
    primer_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc)  {

    if (length(primers) == 0 || nrow(primers) == 0) {
        return(NULL)
    }
    candidates <- NULL  # idx of all selection candidates
    bad.primer.df <- data.frame(Identifier = numeric(0), Index = numeric(0), melting_temp = numeric(0), 
        deltaG_dimer = numeric(0))
    for (idx in seq_len(nrow(primers))) {
        # check whether constraints are fulfilled by selection candidates best primers
        # already found, check if cur primer cvg is the same as the best primer cvg of
        # the current primer is not as high as the selected ones -> don't consider for
        # inclusion
        if (length(candidates) != 0) {
            if (primers[idx, "primer_coverage"] < primers[candidates[1], "primer_coverage"]) {
                break
            }
        }
        if (!is.na(deltaG.cutoff) && length(result) != 0 && nrow(result) != 0) {
            # check for cross-dimerization between current primers and candidate primer
            tmp.result <- suppressWarnings(my_rbind(result, primers[idx, ]))
            comp <- compute.all.cross.dimers(tmp.result, primer_conc, na_salt_conc, 
                    mg_salt_conc, k_salt_conc, tris_salt_conc, 
                    target.temp, check.idx = nrow(tmp.result), 
                    no.structures = TRUE)
            if (any(comp$DeltaG < deltaG.cutoff)) {
                # any deltaG exceeds the allowed cutoff -> don't select this primer
                min.deltaG <- min(comp$DeltaG)
                #message(paste("bad primer id: ", primers[idx, "Identifier"], ", dimerization deltaG: ", 
                  #min.deltaG))
                bad.primer <- data.frame(Identifier = primers[idx, "Identifier"], 
                  Index = idx, melting_temp = NA, deltaG_dimer = min.deltaG)
                bad.primer.df <- rbind(bad.primer.df, bad.primer)
            } else {
                candidates <- c(candidates, idx)
            }
        } else {
            # if there's no results yet, the first primer is automatically a candidate ..
            candidates <- c(candidates, idx)
        }
    }
    if (length(candidates) != 0) {
        # select the best primer from the candidates (all have the same cvg) by
        # resolution through primer efficiency (if available) select primer by highest
        # primer efficiency
        if ("primer_efficiency" %in% colnames(primers)) {
            o <- order(primers[candidates, "mean_primer_efficiency"], decreasing = TRUE)
            idx <- candidates[o][1]
        } else {
            idx <- candidates[1]
        }
    }
    res <- list(sel_idx = idx, filtered_info = bad.primer.df)
    return(res)
}
#' Filter by Melting Temperature Difference
#'
#' Filters primers by melting temperature differences.
#'
#' @param target.temp Target melting temperature in Celsius.
#' @param selected.primers Current candidate primer data frame.
#' @param max.Tm.delta Maximum allowed difference of primer melting temperatures to target temperature.
#' @return Filtered primer data frame.
#' @keywords internal
filter_primers.by.Tm.delta <- function(target.temp, selected.primers, max.Tm.delta) {
    # target.temp: target melting temperature of primers in the optimization set
    # selected.primers: current candidates for selection max.Tm.delta: maximal
    # allowed difference of primer melting temperature to target temperature
    r <- selected.primers
    sel.idx <- NULL
    filtered.df <- data.frame(Identifier = numeric(0), Index = numeric(0), melting_temp = numeric(0), 
        deltaG_dimer = numeric(0))  # excluded primers
    if (!is.na(max.Tm.delta)) {
        Tm.diff <- abs(selected.primers$melting_temp - target.temp)
        sel.idx <- which(Tm.diff > max.Tm.delta)
        if (length(sel.idx) != 0) {
            f.df <- data.frame(Identifier = selected.primers$Identifier[sel.idx], 
                Index = sel.idx, melting_temp = Tm.diff[sel.idx], deltaG_dimer = NA)
            filtered.df <- rbind(filtered.df, f.df)
            r <- selected.primers[-sel.idx, ]
        }
    }
    out <- list(data = r, filtered = Primers(filtered.df))  # data: selected primers, filtered: excluded primers. melting temp entry of filtered: refers to DeltaTm!
    return(out)
}
#' Check of Primer and Template Correspondence.
#'
#' Checks whether the primers relate to the correct templates.
#'
#' @param primer.df An object of class \code{Primers}.
#' @param template.df An object of class \code{Templates}.
#' @return \code{TRUE} if the primers and templates seem to correspond,
#' \code{FALSE} otherwise.
#' @keywords internal
check_correspondence <- function(primer.df, template.df) {
    # check whether primers and templates correspond to one another
    cvd.seqs <- unlist(sapply(strsplit(primer.df$Covered_Seqs, split = ","), function(x) as.numeric(x)))
    if (length(cvd.seqs) == 0) {
        # no covered seqs -> always corresponds
        return(TRUE)
    }
    cvd.seqs <- cvd.seqs[!is.na(cvd.seqs)] # exclude primers for which cvg is not available.
    m <- match(cvd.seqs, template.df$Identifier)
    if (any(is.na(m))) {
        return(FALSE)
        msg <- paste("Could not match all coverage events to the input template data frame.",
            "Please check whether the primers and templates correspond.")
        stop(msg)
    } else {
        return(TRUE)
    }
}
#' Greedy Optimization
#'
#' Greedy approach for solving the primer set coverage problem.
#'
#' @param primers Primer data frame to be optimized.
#' @param template.df Template data frame.
#' @param mode.directionality Primer direction.
#' @param cur.opti.constraints List with optimization constraint settings.
#' @param target.temp Target annealing temperature of the optimized primer set in Celsius.
#' @param allowed.mismatches The number of mismatches primers are allowed to have with the templates.
#' @param opti.limits List with optimization limits.
#' @param primer_conc Primer concentration.
#' @param template_conc Template concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param updateProgress Shiny progress callback function.
#' @return List with optimization data.
#' @keywords internal
optimize.primer.cvg <- function(primers, template.df, mode.directionality, 
    cur.opti.constraints, target.temp, allowed.mismatches, opti.limits, primer_conc, 
    na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, updateProgress = NULL) {

    #print("optimizing primer cvg")
    deltaG.cutoff <- NA
    if ("cross_dimerization" %in% names(cur.opti.constraints)) {
        deltaG.cutoff <- cur.opti.constraints[["cross_dimerization"]]["min"]
    }
    l.df <- template.df
    max.cvg <- get_cvg_ratio(primers, template.df)
    cur.cvg <- 0
    it <- 0
    result <- NULL # not initialized here due to Primers class.
    all.bad.primers <- data.frame(Identifier = numeric(0), Index = numeric(0), 
                        melting_temp = numeric(0), deltaG_dimer = numeric(0))
    # order primers by cvg for greedy selection
    selected.primers <- primers[order(primers$primer_coverage, decreasing = TRUE), ]  
    while ((cur.cvg != max.cvg) && nrow(selected.primers) >= 1) {
        it <- it + 1
        #message(it)
        best.primer.data <- select.best.primer.idx(result, selected.primers, 
            deltaG.cutoff, target.temp, primer_conc, 
            na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc)
        if (length(best.primer.data) == 0 || is.na(best.primer.data$sel_idx)) {
            # no primer could be selected..
            selected.primers <- primers[FALSE, ]
            break
        } else {
            best.primer.idx <- best.primer.data$sel_idx
            bad.primer.df <- best.primer.data$filtered_info  # e.g. cross-dimerizing primers
            all.bad.primers <- rbind(all.bad.primers, bad.primer.df)
            best.covered <- as.numeric(unlist(strsplit(selected.primers[best.primer.idx, "Covered_Seqs"], split = ",")))
            if (length(best.covered) == 0) {
                # cannot improve anymore: best (highest cvg primer) didn't cover anything
                break
            }
            # selected.primers gives only the coverage count for sequences that haven't been covered yet :-)
            cumulative.coverage <- sum(result$Marginal_Coverage_Gain) + selected.primers[best.primer.idx, "primer_coverage"]  
            # compute current coverage
            cur.cvg <- cumulative.coverage/nrow(template.df)
            # use primer entry from original table to have the 'true' covered seqs in the
            # result, not only the marginal gain (primer_coverage is not the marginal gain!!)
            primer.entry <- primers[match(selected.primers[best.primer.idx, "Identifier"], primers$Identifier), ]
            entry <- data.frame(Iteration = it, 
                                Marginal_Coverage_Gain = selected.primers[best.primer.idx, "primer_coverage"], 
                                Cumulative_Coverage = cumulative.coverage, Cumulative_Coverage_Ratio = cur.cvg, 
                                primer.entry, stringsAsFactors = FALSE)
            result <- rbind(result, entry)
            rm.idx <- c(bad.primer.df$Index, best.primer.idx)
            selected.primers <- selected.primers[-rm.idx, ]  # remove the selected primer and selection candidates that were filtered this round
            # re-compute coverage by removing templates that are covered by newly added primer (best.covered)
            selected.primers <- evaluate.diff.primer.cvg(selected.primers, best.covered, template.df)  
        }
    }
    # reset cvg to consistent values from before (not marginal gain anymore):
    if (length(result) != 0) {
        m <- match(result$Identifier, primers$Identifier)
        result[, colnames(primers)] <- primers[m, ]
        # remove redundant primers
        cvg.matrix <- get.coverage.matrix(primers, template.df)
        primer.idx <- remove.redundant.cols(m, cvg.matrix)  # updated primer.idx with unnecessary primers removed
        result <- result[match(primers$Identifier[primer.idx], result$Identifier), , drop = FALSE]
    } else {
        result <- primers # could not optimize anything
    }
    out <- list(data = Primers(result), bad_primers = all.bad.primers)
    return(out)
}
#' Greedy Optimization.
#'
#' Greedy approach for solving the primer set coverage problem.
#'
#' @param primers Primer data frame to be optimized.
#' @param settings A \code{DesignSettings} object.
#' @param template.df Template data frame.
#' @param mode.directionality Primer direction.
#' @param required.cvg Target coverage ratio.
#' @param allowed.mismatches The number of mismatches primers are allowed to have with the templates.
#' @param primer_conc Primer concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param template_conc Template data frame.
#' @param allowed.other.binding.ratio Ratio of primers allowed to bind to non-target regions.
#' @param allowed.stop.codons  Consider mismatch binding events that induce stop codons.
#' @param allowed.region.definition Definition of the target binding sites used for evaluating the coverage.
#' If \code{allowed.region.definition} is \code{within}, primers have to lie within the allowed binding region.
#' If \code{allowed.region.definition} is \code{any}, primers have to overlap with the allowed binding region.
#' The default is that primers have to bind within the target binding region.
#' @param disallowed.mismatch.pos The number of positions from the primer 3' end where mismatches should not be allowed.
#' All primers binding templates with mismatches within \code{disallowed.mismatch.pos} from the 3' end are disregarded.
#' @param target.temps Target melting temperatures for optimized sets in Celsius.
#' @param fw.primers List with already optimized primer data frames corresponding to \code{target.temps}.
#' @param updateProgress Shiny progress callback function.
#' @return List with optimization data.
#' @keywords internal
select.primers.by.cvg <- function(primers, settings, template.df, mode.directionality = c("fw", "rev"), 
    required.cvg = 1, allowed.mismatches, primer_conc, na_salt_conc, mg_salt_conc, 
    k_salt_conc, tris_salt_conc, template_conc, allowed.other.binding.ratio,
    allowed.stop.codons, allowed.region.definition = c("within", "any"), 
    disallowed.mismatch.pos, target.temps = NULL, fw.primers = NULL, updateProgress = NULL) {

    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' arg.")
    }
    mode.directionality <- match.arg(mode.directionality)
    if (length(allowed.region.definition) == 0) {
        stop("Please supply the 'allowed.region.definition' arg.")
    }
    allowed.region.definition <- match.arg(allowed.region.definition)
    opti.constraints <- opti(settings)
    opti.limits <- optiLimits(settings)
    if (length(primers) == 0) {
        return(NULL)
    } else if (nrow(primers) == 0) {
        # nothing to optimize
        return(list(opti = primers, all_results = list(primers), all_used_constraints = list(settings), 
            used_constraints = settings))
    }
    # add primer cvg if not present yet: needed for optimization
    if (!"primer_coverage" %in% colnames(primers)) {
        stop("Optimization requires pre-computed primer coverage")
    }
    # determine computation constraints:
    Tm.brackets <- create.Tm.brackets(primers, template.df, settings, target.temps)
    primers <- Tm.brackets$primers
    Tm.groups <- Tm.brackets$group_idx
    Tm.group.df <- Tm.brackets$df
    opti.results <- vector("list", nrow(Tm.group.df))
    # 1. precompute constraints for all temperatures
    diagnostic.location <- NULL
    Tm.data <- compute.Tm.sets(primers, template.df, Tm.brackets, settings, mode.directionality, 
        primer_conc, template_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, 
        allowed.mismatches, allowed.other.binding.ratio, 
        allowed.stop.codons, allowed.region.definition, disallowed.mismatch.pos, 
        TRUE, required.cvg,
        fw.primers, diagnostic.location, updateProgress)
    Tm.sets <- Tm.data$sets
    Tm.settings <- Tm.data$settings
    #############
    # 2. optimize all temperature primer sets
    Tm.group.idx <- NULL
    max.nbr.relaxations <- 3
    opti.results.data <- vector("list", length(Tm.sets)) # for debug
    for (Tm.group.idx in seq_along(Tm.sets)) { # for debug only
    #opti.results.data <- foreach(Tm.group.idx = seq_along(Tm.sets), .combine = c) %dopar% {
        Tm.set <- Tm.sets[[Tm.group.idx]]
        cur.settings <- Tm.settings[[Tm.group.idx]]
        if (length(Tm.set) == 0 || nrow(Tm.set) == 0) {
            # no Tm set -> nothing to optimize. this shouldn't really happen except if there's no data.
            Tm.out <- list(list(data = Tm.set, used_constraints = opti.constraints))
        } else {
            #message(paste("Nbr of primers: ", dim(Tm.set)[1], sep = ""))
            target.temp <- Tm.group.df$target_Tm[Tm.group.idx]
            annealing.temp <- Tm.group.df$annealing_temp[Tm.group.idx]
            if (is.function(updateProgress)) {
                detail <- paste("Tm group: ", Tm.group.idx, sep = "")
                updateProgress(1/length(opti.results), detail, "inc")
            }
            #message(paste("Optimizing for Tm-group: ", Tm.group.idx, "/", nrow(Tm.group.df), 
            #    sep = ""))
            target.cvg <- min(get_cvg_ratio(Tm.set, template.df), required.cvg)
            result.data <- optimize.primer.cvg(Tm.set, template.df, mode.directionality, 
                opti(cur.settings), annealing.temp, allowed.mismatches, 
                opti.limits, primer_conc, na_salt_conc, mg_salt_conc, k_salt_conc, 
                tris_salt_conc, updateProgress)
            result <- result.data$data
            result.cvg <- get_cvg_ratio(result, template.df)
            # relaxation
            while (result.cvg < target.cvg && "cross_dimerization" %in% names(opti(cur.settings))) {
                # don't lose any coverage here in the optimization ...
                # adjust DeltaG
                message("Greedy deltaG relaxation. Current coverage: ", result.cvg, 
                            "Target coverage: ", target.cvg, 
                            ". Old deltaG: ", constraints(cur.settings)$cross_dimerization["min"])
                constraints(cur.settings)$cross_dimerization <- relax.opti.constraints(opti(cur.settings), optiLimits(Tm.settings[[Tm.group.idx]]), 
                                                                                      opti(Tm.settings[[Tm.group.idx]]))$cross_dimerization
                constraintLimits(cur.settings)$cross_dimerization <- set.new.limits(optiLimits(cur.settings), optiLimits(Tm.settings[[Tm.group.idx]]), 
                                                                                      opti(Tm.settings[[Tm.group.idx]]))$cross_dimerization
                result.data <- optimize.primer.cvg(Tm.set, template.df, mode.directionality, 
                    opti(cur.settings), annealing.temp, allowed.mismatches, 
                    opti.limits, primer_conc, na_salt_conc, mg_salt_conc, k_salt_conc, 
                    tris_salt_conc,  updateProgress)
               result <- result.data$data
               result.cvg <- get_cvg_ratio(result, template.df)
               message("-> New coverage: ", result.cvg, 
                            ". New deltaG: ", constraints(cur.settings)$cross_dimerization["min"])
            }
            Tm.out <- list(list(data = result.data$data, used_constraints = cur.settings))
            opti.results.data[[Tm.group.idx]] <- Tm.out[[1]] # for debug only
        }
    }
    opti.results <- lapply(opti.results.data, function(x) x$data)  # extract optimized sets
    names(opti.results) <- Tm.group.df$target_Tm # uses the target Tm's before relaxation
    best.idx <- select.best.opti.result(opti.results, template.df)
    if (length(best.idx) == 0) {
        # could not determine a solution, arbitrarily return the first opti.result
        best.idx <- 1
    }
    opti.result <- opti.results[[best.idx]]
    all.used.constraints <- lapply(opti.results.data, function(x) x$used_constraints)  # all used settings
    result <- list(opti = opti.result, 
        used_constraints = opti.results.data[[best.idx]]$used_constraints, 
        all_used_constraints = all.used.constraints, 
        all_results = opti.results)
    #message("optimization w/ greedy finished")
    return(result)
}
#' Selection of Best Greedy Result
#'
#' Selects best greedy primer data set.
#'
#' @param opti.results List with primer data frames for different target melting temperatures.
#' @param template.df Template data frame.
#' @return The index of the best primer data set.
#' @keywords internal
select.best.opti.result <- function(opti.results, template.df) {
    # could have a single function for ILP and greedy result selection?
    cvg.per.group <- sapply(opti.results, function(x) if (length(x) != 0 && nrow(x) != 
        0) {
        max(x$Cumulative_Coverage)
    } else {
        0
    })
    if (all(cvg.per.group == 0)) {
        message("Could not find a solution with primers covering any of the templates.")
        if (length(opti.results) > 0) {
            return(1)
        } else {
            return(NULL)
        }
    }
    cvg.ratios <- sapply(cvg.per.group, function(x) x/nrow(template.df))
    #message("Cumulative coverage per group: ")
    #message(cvg.per.group)
    #message("Coverage ratios: ")
    #message(cvg.ratios)
    # TIE BREAKING: at the moment: just use primer_efficiency: later: score primer
    # set using all available parameters
    eps <- 0.01  # if the cvg difference is just 1%, choose the set with better efficiency
    m.cvg <- max(cvg.per.group, na.rm = TRUE)
    max.idx <- which(cvg.per.group == m.cvg)[1]
    sel.idx <- which((cvg.ratios[max.idx] - cvg.ratios) < eps)
    idx <- NULL
    if ("primer_efficiency" %in% colnames(opti.results[[1]])) {
        eff <- unlist(lapply(opti.results[sel.idx], function(x) mean(x$mean_primer_efficiency, 
            na.rm = TRUE)))
        # message(opti.results[[1]]$primer_efficiency)
        #message(paste("Primer set mean efficiencies:", eff, sep = ""))
        idx <- sel.idx[which.max(eff)]
    } else {
        # no primer efficiency available -> just select the max
        idx <- max.idx[1]
    }
    #message(paste("Selected primer set: ", idx, sep = ""))
    return(idx)
}


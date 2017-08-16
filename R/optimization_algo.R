############
# Primer optimization algorithms
###########

#' File Name for Initialized Primers.
#'
#' Constructs a filename for initialized primers.
#'
#' @param cur.results.loc Directory where the file should be stored.
#' @param GROUP Sample name of templates.
#' @param primer.lengths Interval of desired primer lengths.
#' @param mode.directionality Directionality of the primers
#' @param allowed.region.definition Definition of the allowed region.
#' @param init.algo Initialization algorithm identifier.
#' @param max.degen Maximum degeneracy of primers.
#' @param conservation Required ratio of primer conservation.
#' @return A filename for the initialized primers.
#' @keywords internal
get.init.file.name <- function(cur.results.loc, GROUP, primer.lengths, mode.directionality, 
    allowed.region.definition, init.algo, max.degen, conservation) {
    # get the filename for the initial primer set
    if (length(cur.results.loc) == 0) {
        return(NULL)
    }
    algo.info <- ""
    if (init.algo == "tree") {
        algo.info <- paste(algo.info, "degen=", max.degen, "_con=", conservation, 
            sep = "")
    }
    path <- file.path(cur.results.loc, paste(GROUP, "_initial_primers_", init.algo, 
        "_", primer.lengths[1], "_", primer.lengths[2], "_", algo.info, mode.directionality, 
        ".csv", sep = ""))
    return(path)
}

#' Write Out Optimization Data
#'
#' Writes out all data relating to the optimization of primers.
#'
#' @param opti.results.loc Folder where optimization data reside.
#' @param optimal.primers.data List with optimization results.
#' @param mode.directionality Direction of primers.
#' @param settings Settings used in the optimization procedure.
#' List containing fw, rev settings.
#' @param sample.name Name of template sample.
#' @param template.df Template data frame.
#' @param max.degen Maximal degeneracy of primers.
#' @return Write-out of primer information to \code{opti.results.loc}.
#' @keywords internal
write.out.primer.info <- function(opti.results.loc, optimal.primers.data, 
    mode.directionality, settings, sample.name, template.df, max.degen) {
    if (length(opti.results.loc) != 0) {
        optimal.primers <- optimal.primers.data$opti
        all.opti.results <- optimal.primers.data$all_results
        for (i in seq_along(all.opti.results)) {
            write.csv(all.opti.results[[i]], file = file.path(opti.results.loc, paste(sample.name, 
                "_optimized_set_target_temp_", names(all.opti.results)[i], ".csv", 
                sep = "")), row.names = FALSE)
        }
        template.df <- update_template_cvg(template.df, optimal.primers)
        write.csv(template.df, file = file.path(opti.results.loc, paste(sample.name, "_optimized_templates.csv", 
            sep = "")), row.names = FALSE)
        # save constraints:
        for (i in seq_along(settings)) {
            setting <- settings[[i]]
            write_settings(setting, file.path(opti.results.loc, paste0(sample.name, 
                           "_constraint_settings_", names(settings)[i], ".xml")))
        }
        write.csv(optimal.primers, file = file.path(opti.results.loc, paste(sample.name, 
            "_FINAL_optimized_primers.csv", sep = "")), row.names = FALSE)
        merged.result <- merge.ambig.primers(optimal.primers, mode.directionality, 
            max.degen)
        write_primers(merged.result, file.path(opti.results.loc, 
                        paste(sample.name, "_FINAL_optimized_primers.fasta", sep = "")))
    }
}
#' Design of Multiplex PCR Primers.
#'
#' Designs a primer set maximizing the number of covered templates using
#' the smallest possible number of primers. The algorithm tries to ensure
#' that the designed set of primers achieves a coverage ratio not lower than
#' \code{required.cvg}. To this end, the constraints for designing
#' primers may be relaxed.
#'
#' @section{1. Initialization}:
#' The primer design algorithm consists
#' of three steps: primer initialization, filtering, and optimization.
#' The method for initializing a set of candidate primers is determined
#' via \code{init.algo}. If \code{init.algo} is set to \emph{naive}, primers
#' are created by extracting substrings from all input template sequences.
#' If \code{init.algo} is set to \emph{tree}, degenerate primers are created by
#' merging similar subsequences by forming their consensus sequence up to
#' a degeneracy of at most \code{max.degen}. The tree-based initialization
#' is recommended for related sequences.
#'
#' @section{2. Filtering}:
#' The candidate primer set is filtered according to the constraints
#' specified in the \code{settings} object. In some cases, it is necessary
#' to relax the constraints in order to reach the desired \code{required.cvg}.
#' In these cases, primers that fail the input constraints may be selected. 
#' If you would like to skip the initialization and filtering stages,
#' you can provide an evaluated \code{Primers} object via \code{primer.df}.
#'
#' @section{3. Optimization}:
#' Optimizing a primer set entails finding the smallest subset of primers
#' maximizing the coverage, which is done by solving the set cover problem.
#' If melting temperature differences are a constraint,
#' the optimization procedure automatically samples ranges of melting
#' temperatures to find optimal sets for all possible temperatures.
#' You can select the used optimization algorithm via \code{optia.algo}, where
#' you can set "Greedy" for a greedy algorithm or "ILP for 
#' an integer linear program formulation (ILP).
#' While the worst-case runtime of the
#' greedy algorithm is shorter than the worst-case runtime of the ILP, 
#' the greedy solution may yield larger primer sets than the ILP solution.
#'
#' @param template.df A \code{Templates} object containing the template
#' sequences and target regions for designing primers.
#' @param mode.directionality The template strand for which primers shall be designed.
#' Primers can be designed either for forward strands ("fw"), 
#' for reverse strands ("rev"), or for both strands ("both"). The default setting
#' is "both".
#' @param settings A \code{DesignSettings} object specifying the constraint settings for filtering and optimization.
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
#' Only applicable for the tree-based primer initialization. At the its
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
#' @return A list with the following fields:
#' \describe{
#' \item{\code{opti}:}{A \code{Primers} object providing the designed primer set.}
#' \item{\code{used_constraints}:}{A list with \code{DesignSettings} objects
#' for each primer direction providing the (possibly relaxed) constraints used
#' for designing the optimal primers.}
#' \item{\code{all_results}:}{A list containing objects of class \code{Primers}.
#' Each list entry corresponds to an optimal primer set for a given
#' melting temperature.}
#' \item{\code{all_used_constraints}:}{A list containing \code{DesignSettings} object for each optimized set in \code{all_results}.}
#' \item{\code{filtered}:}{A list with data providing information on the results
#' of the filtering procedure.}
#' } 
#' @family primer functions
#' @note Some constraints specified in the \code{settings} object
#' can only be computed if additional software is installed,
#' please see the documentation of \code{\link{DesignSettings}} for an overview of all possible settings and the \code{\link{ConstraintSettings}} documentation
#' for an overview of all possible constraints.
#' Usage of \code{init.algo = "tree"} requires an installation of
#' the multiple alignment program MAFFT (http://mafft.cbrc.jp/alignment/software/).
#' @export
#' @keywords Primers
#' @examples
#' # Define PCR settings and primer criteria
#' data(Ippolito)
#' constraints(settings)$primer_length <- c("min" = 18, "max" = 18)
#' # Design only forward primers using a greedy algorithm
#' optimal.primers.greedy <- design_primers(template.df[1:2,], "fw", settings, init.algo = "naive")
#' # Design forward and reverse primers using an ILP and store 
#' # the results in 'out.dir'
#' out.dir <- tempdir()
#' optimal.primers.ILP <- design_primers(template.df[1:2,], "both", settings,
#'                          init.algo = "naive", opti.algo = "ILP",
#'                          cur.results.loc = out.dir)
#' # Usage of the tree-based initialization strategy (requires MAFFT)
#' \dontrun{
#' optimal.primers.tree <- design_primers(template.df[1:2,], "both", settings,
#'                          init.algo = "tree", opti.algo = "ILP",
#'                          max.degen = 16,
#'                          cur.results.loc = out.dir)
#' }
design_primers <- function(template.df, mode.directionality = c("both", "fw", "rev"),
     settings, init.algo =  c("naive", "tree"), opti.algo = c("Greedy", "ILP"),
     required.cvg = 1.0, timeout = Inf, max.degen = 16, conservation = 1.0,
     sample.name = NULL, cur.results.loc = NULL, 
     primer.df = NULL, updateProgress = NULL) {
  
    # never allow any binding outside the target region for designing
    if (init.algo == "tree" && !check.tool.function()["MAFFT"]) {
        stop("If you would like to use the tree-based primer initialization strategy, please install MAFFT (http://mafft.cbrc.jp/alignment/software/) first.")
    }

    conOptions(settings)$allowed_other_binding_ratio <- 0
    mode.directionality <- match.arg(mode.directionality)
    opti.algo <- match.arg(opti.algo)
    init.algo <- match.arg(init.algo)
    if (length(settings) == 0 || !is(settings, "DesignSettings")) {
        stop("Please provide a DesignSettings object.")
    }
    if (!is(template.df, "Templates")) {
        stop("Please supply a valid template data frame.")
    }
    if (required.cvg < 0 || required.cvg > 1) {
        stop("Required coverage should be in [0,1].")
    }
    if (max.degen < 0) {
        stop("The maximal degeneracy should be positive.")
    }
    if (conservation < 0 || conservation > 1) {
        stop("The top-conservation percentile should be in [0,1].")
    }
    # check whether all required constraints are provided
    required.constraints <- c("primer_length", "primer_coverage")
    if (!all(required.constraints %in% names(constraints(settings)))) {
        stop("Please provide a settings object containing at least primer_length and primer_coverage as constraints.")
    }
    if (is.null(sample.name)) {
        if (nrow(template.df) >= 1)  {
            sample.name <- template.df$Run[1]
        }
    }
    if (sample.name == "") {
        # sample name should be different from "" for input from csv
        sample.name <- "openPrimeR"
    }
    allowed.region.definition <- conOptions(settings)$allowed_region_definition
    # check whether binding regions are long enough to design primers
    length.check <- check.init.primer.length(template.df, allowed.region.definition,
                    filters(settings)$primer_length, mode.directionality)
    if (!length.check) {
        # allowed region is too short
        my.warning("AllowedRegionTooShort", "The allowed binding region was too short. Please adjust the binding region such that primers can be created for all templates.")
        return(NULL)
    }
    if (length(primer.df) != 0) {
        # skip the init/filtering steps
        if (!is(primer.df, "Primers")) {
            stop("primer.df is not a Primers objectc.")
        } 
        if (!"primer_coverage" %in% colnames(primer.df)) {
            stop("primer.df does not have coverage entries.")
        }
        fw.primers <- primer.df[primer.df$Forward != "", ]
        rev.primers <- primer.df[primer.df$Reverse != "", ]
        if (mode.directionality == "fw") {
            single.primers <- fw.primers
        } else if (mode.directionality == "rev") {
            single.primers <- rev.primers
        } else {
            if (nrow(fw.primers) == 0 || nrow(rev.primers) == 0) {
                warning("'mode.directionality' is 'both', but only primers of one directionality were provided.") 
            }
            single.primers <- NULL
        }
    } else {
        fw.primers <- NULL
        rev.primers <- NULL
        single.primers <- NULL
    }
    target.temps <- NULL
    if (mode.directionality == "fw" || mode.directionality == "rev") {
        optimal.primers.data <- design_primers.single(template.df, sample.name,
            mode.directionality, settings, timeout, opti.algo, 
            allowed.region.definition, init.algo, max.degen, 
            conservation, target.temps, required.cvg, 
            cur.results.loc = cur.results.loc, 
            primer.df = single.primers, updateProgress = updateProgress)
    } else {
        # both directions are to be optimized
        message("#####\n# (FW) Primer design for forward primers\n#####")
        optimal.primer.data.fw <- design_primers.single(template.df, sample.name, "fw", 
            settings, timeout, opti.algo, allowed.region.definition, init.algo, max.degen, 
            conservation, target.temps, required.cvg, 
            cur.results.loc = cur.results.loc,
            primer.df = fw.primers, updateProgress = updateProgress)
        if (is.function(updateProgress)) { # reset progress to 0 for rev opti
            updateProgress(0, "", "set")
        }
        # exclude primer temperature-sets with 0 cvg
        non.zero.cvg <- unlist(lapply(optimal.primer.data.fw$all_results, function(x) sum(x$primer_coverage) != 0))
        optimal.primer.data.fw$all_results <- optimal.primer.data.fw$all_results[non.zero.cvg]
        optimal.primer.data.fw$all_used_constraints <- optimal.primer.data.fw$all_used_constraints[non.zero.cvg]
        if (length(optimal.primer.data.fw$all_results) == 0) {
            warning("Could not construct a primer set for both directions with any coverage given the current settings.")
            return(NULL)
        }
        # optimize reverse primers for the same target temperatures as the fw primers
        if ("melting_temp_diff" %in% names(opti(settings))) {
            target.temps <- unlist(sapply(seq_along(optimal.primer.data.fw$all_results), function(x) mean(optimal.primer.data.fw$all_results[[x]][,"melting_temp"], na.rm= TRUE)))
            names(optimal.primer.data.fw$all_results) <- target.temps
            if ("melting_temp_range" %in% names(filters(settings))) {
                # adjust the melting_temp range here to save time later
                delta <- unname(opti(settings)$melting_temp_diff["max"] / 2)
                new.Tm.range <- c("min" = min(target.temps) - delta, "max" = max(target.temps) + delta)
                new.Tm.limit <- c(new.Tm.range["min"] - delta, new.Tm.range["max"] + delta)
                constraints(settings)$melting_temp_range <- new.Tm.range
                constraintLimits(settings)$melting_temp_range <- new.Tm.limit
            }
        } else {
            target.temps <- rep(NA, length(optimal.primer.data.fw$all_results))
        }
        # adjust the templates for rever design to include only templates that can be covered by both direction primers
        opti.fw <- optimal.primer.data.fw$opti
        cvg.idx <- unique(unlist(covered.seqs.to.idx(opti.fw$Covered_Seqs, template.df)))
        required.nbr <- required.cvg * nrow(template.df) # required nbr of covered templates
        rev.required.cvg <- min(length(cvg.idx) / required.nbr, 1) # update required cvg to the new value considering only the already-covered templates
        rev.template.df <- template.df[cvg.idx[order(cvg.idx)],] # update template data frame to retain only the already covered templates
        primers.fw <- optimal.primer.data.fw$all_results  # arg for cross dimerization considerations between fw/rev primers
        #print(target.temps)
		message("#####\n# (REV) Primer design for reverse primers\n#####")
        optimal.primer.data.rev <- design_primers.single(rev.template.df, sample.name, "rev", 
            settings, timeout, opti.algo, allowed.region.definition, init.algo,
            max.degen, conservation, target.temps, rev.required.cvg, 
            cur.results.loc = cur.results.loc, primer.df = rev.primers, updateProgress = updateProgress)
		message("#####\n# (BOTH) Aggregating results\n#####")
        # TODO: read.table error here after designing for 'both' in test program
        opti.fw <- optimal.primer.data.fw$all_results # base consideration of templates on the best set from the 'fw' run
        opti.rev <- optimal.primer.data.rev$all_results
        #print("FW set:")
        #print(opti.fw)
        #print("REV set:")
        #print(opti.rev)
        if ("melting_temp_diff" %in% names(constraints(settings))) {
            # match fw and rev primers for melting temps if melting temp diff is active
            # matching of sets should be improved ...?
            allowed.diff.fw <- max(sapply(optimal.primer.data.fw$all_used_constraints, function(x) constraints(x)$melting_temp_diff))
            allowed.diff.rev <- max(sapply(optimal.primer.data.fw$all_used_constraints, function(x) constraints(x)$melting_temp_diff))
            allowed.diff <- max(c(allowed.diff.fw, allowed.diff.rev))
            fw.tm <- as.numeric(names(opti.fw))
            rev.tm <- unlist(sapply(seq_along(optimal.primer.data.rev$all_results), function(x) mean(optimal.primer.data.rev$all_results[[x]][,"melting_temp"], na.rm= TRUE)))
            compatible <- do.call(rbind, lapply(fw.tm, function(x) abs(x - rev.tm) < allowed.diff))
            opti.fw.indices <- unlist(lapply(seq_len(nrow(compatible)), function(x)
                rep(x, length(which(compatible[x,])))))
            opti.fw <- opti.fw[opti.fw.indices]
            opti.rev.indices <- unlist(lapply(seq_len(ncol(compatible)), function(x)
                rep(x, length(which(compatible[,x])))))
            opti.rev <- opti.rev[opti.rev.indices]
        }
        fw.rev.data <- evaluate.fw.rev.combinations(opti.fw, opti.rev, template.df, target.temps)
        sel.set <- select.best.primer.set(fw.rev.data$stats)
        if (length(sel.set) == 0) {
            warning("Could not select an optimal fw-rev combination of primers.")
            optimal.primers <- NULL
        } else {
            optimal.primers <- fw.rev.data$sets[[sel.set]]
        }
        # update melting temp diff and cross dimerization to account for fw/rev primers
        for (i in seq_along(fw.rev.data$sets)) {
            cur.set <- fw.rev.data$sets[[i]]
            # TODO: update.opti.results -> Error in read.table
            cur.set <- update.opti.results(cur.set, settings, template.df)
            fw.rev.data$sets[[i]] <- cur.set
        }
        cur.set <- optimal.primers
        cur.set <- update.opti.results(cur.set, settings, template.df)
        optimal.primers <- cur.set
        # construct result: settings
        if (length(sel.set) != 0) {
            used.settings.fw <- optimal.primer.data.fw$all_used_constraints[[sel.set]]
            used.settings.rev <- optimal.primer.data.rev$all_used_constraints[[sel.set]]
        } else {
            used.settings.fw <- NULL
            used.settings.rev <- NULL
        }
        used.settings <- list("fw" = used.settings.fw, "rev" = used.settings.rev)
        # unify filtering stats from fw and rev primers
        filter.data <- list()
        for (i in seq_along(optimal.primer.data.fw$filtered)) {
            field <- names(optimal.primer.data.fw$filtered)[i]
            if (field != "used_settings") {
                # constraints are dealt with in another way ..
                filter.data[[field]] <- my_rbind(optimal.primer.data.fw$filtered[[i]], optimal.primer.data.rev$filtered[[i]])
            }
        }
        # gather filtered primers: combine fw and rev filtered primers
        optimal.primers.data <- list(opti = optimal.primers, all_results = fw.rev.data$sets, 
            used_constraints = used.settings,
            filtered = filter.data)
        if (length(cur.results.loc) != 0 && length(sel.set) != 0) {
            # store data
            message("Storing results")
            cur.results.loc <- file.path(cur.results.loc, "both")
            dir.create(cur.results.loc, showWarnings = FALSE)
            PCR.settings <- PCR(settings)
            template.df <- update_template_cvg(template.df, optimal.primers)
            write.out.primer.info(cur.results.loc, optimal.primers.data, 
                mode.directionality, used.settings,
                sample.name, template.df, max.degen)
            plot.all.cvg.info(sample.name, cur.results.loc, optimal.primers, template.df, 
                mode.directionality, "optimized", PCR.settings$primer_concentration,
                PCR.settings$Na_concentration, PCR.settings$Mg_concentration, 
                PCR.settings$K_concentration, PCR.settings$Tris_concentration, 
                settings, required.cvg, used.settings)  # plot opti coverage info
        }
    }
    opti.data <- optimal.primers.data$all_results
    cvg <- unlist(lapply(opti.data, function(x) get_cvg_ratio(x, template.df, mode.directionality = mode.directionality, as.char = TRUE)))
    message("Coverage of each designed set: ", paste(cvg, collapse = ", "))
    return(optimal.primers.data)
}
#' Selection of Best Primer Set.
#' 
#' Selects the best primer set according to coverage and melting temperature differences among primers in the set.
#'
#' @param stats Statistics of the primer sets to be evaluated.
#' @return The index of the best primer set.
#' @keywords internal
select.best.primer.set <- function(stats) {
    if (length(stats) == 0) {
        return(NULL)
    }
    # stats: properties of each considered primer set
    idx <- which(stats$Coverage == max(stats$Coverage))  # highest coverage sets
    if (any(!is.na(stats$TempDiff[idx]))) {
        idx <- idx[which.min(stats$TempDiff[idx])]  # highest cvg set with lowest temperature difference between primers
    }
    return(idx)
}
#' Evaluation of Set Combinations
#'
#' Evaluates the combinations of forward and reverse primer sets.
#'
#' @param opti.fw List with optimized forward primer sets.
#' @param opti.rev List with optimized reverse primer sets.
#' @param template.df Template data frame for which primers were designed.
#' @param cur.target.temps Target melting temperatures for each primer data frame given
#' in \code{opti.fw} and \code{opti.rev}.
#' @return List with information on the combinations of forward and reverse primers
#' as well as the combined data frames themselves.
#' @keywords internal
evaluate.fw.rev.combinations <- function(opti.fw, opti.rev, template.df, cur.target.temps) {
    # given optimized fw and rev primers with corresponding target temperatures,
    # determines properties of the corresponding sets
    # opti.fw: list with data frames giving fw optimized primers opti.rev: list with
    # data frames giving rev optimized primers template.df: templates
    if (length(opti.fw) != length(opti.rev)) {
        stop("Forward/reverse primer sets did not correspond to each other.
            Number of fw primer sets: ", length(opti.fw), "\n",
            "Number of rev primer sets: ", length(opti.rev))
    }
    stat.df <- NULL
    primer.sets <- vector("list", length(opti.fw))
    # set fw temperatures as names as temperatures of fw/rev sets should correspond:
    names(primer.sets) <- names(opti.fw)
    for (i in seq_along(opti.fw)) {
        combi.df <- my_rbind(opti.fw[[i]], opti.rev[[i]])  # better than my_rbind
        primer.sets[[i]] <- combi.df
        cvg <- get_cvg_ratio(combi.df, template.df)
        temp.diff <- NA
        if ("melting_temp" %in% colnames(combi.df)) {
            temp.diff <- max(combi.df$melting_temp) - min(combi.df$melting_temp)
        }
        stats <- data.frame(Target_Temperature = cur.target.temps[i], Coverage = cvg, 
            TempDiff = temp.diff)
        stat.df <- rbind(stat.df, stats)
    }
    out <- list(stats = stat.df, sets = primer.sets)
    return(out)
}

#' Augmentation of Optimized Primer Data.
#'
#' Adds melting_temp_diff and cross_dimerization info to optimized sets.
#'
#' @param primer.df A primer data frame.
#' @param settings A \code{DesignSettings} object.
#' @return An updated primer data frame.
#' @keywords internal
update.opti.results <- function(primer.df, settings, template.df) {
    if (length(primer.df) == 0 || nrow(primer.df) == 0) {
        return(primer.df)
    }
    mode.directionality <- get.analysis.mode(primer.df)
    annealing.temp <- compute_annealing_temp(primer.df, mode.directionality, 
                                template.df, PCR(settings)$Na_concentration,
                                PCR(settings)$Mg_concentration, 
                                PCR(settings)$K_concentration, 
                                PCR(settings)$Tris_concentration, 
                                PCR(settings)$primer_concentration)
    if ("melting_temp_diff" %in% colnames(primer.df)) {
        temp.diff <- get.melting.temp.diff(primer.df$Tm_C_fw, primer.df$Tm_C_rev)
        primer.df[, "melting_temp_diff"] <- temp.diff
        primer.df[, "EVAL_melting_temp_diff"] <- temp.diff <= opti(settings)$melting_temp_diff["max"]
    }
    cross.dimer.df <- NULL
    if ("cross_dimerization" %in% names(constraints(settings))) {
        # here we have no_structs not set, since we want to annotate the designed data correctly
        cross.dimer.df <- compute.all.cross.dimers.frontend(primer.df, PCR(settings)$primer_concentration, 
                    PCR(settings)$Na_concentration, PCR(settings)$Mg_concentration, 
                    PCR(settings)$K_concentration, PCR(settings)$Tris_concentration, 
                    annealing.temp)
    }
    if (length(cross.dimer.df) != 0) {
        # no entries in constraint.values when there's nothing to report ...
        primer.df <- update.constraint.values(primer.df, cross.dimer.df)
        primer.df$EVAL_cross_dimerization <- (primer.df$Cross_Dimer_DeltaG >= constraints(settings)$cross_dimerization["min"])
    }
    return(primer.df)
}
#' Design Primers for a Single Direction
#'
#' Designs primers for a single direction.
#'
#'
#' @param template.df Template data frame with sequences for which primers shall be designed.
#' @param sample.name Identifier for the templates.
#' @param mode.directionality Template strands for which primers shall be designed.
#' Primers can be designed either only for forward strands, only for reverse strands, or
#' both strand directions.
#' @param settings A \code{DesignSettings} object specifying the 
#' criteria for designing primers.
#' @param timeout Timeout in seconds for the optimization with ILPs.
#' @param opti.algo The algorithm to be used for solving the primer set covering problem.
#' @param allowed.region.definition Definition of the target binding sites used for evaluating the coverage.
#' If \code{allowed.region.definition} is "within", primers have to lie within the allowed binding region.
#' If \code{allowed.region.definition} is "any", primers have to overlap with the allowed binding region.
#' The default is that primers have to bind within the target binding region.
#' @param init.algo The algorithm to be used for initializing primers.
#' If \code{init.algo} is \code{naive}, then primers are constructed from substrings of the input template sequences.
#' If \code{init.algo} is \code{tree}, phylogenetic trees are used to form degenerate primers whose degeneracy is
#' bounded by \code{max.degen}. 
#' @param max.degen Maximal degeneracy of merged primers.
#' @param conservation When using the tree-based primer initialization, consider only the \code{conservation}
#' percentile of regions with the highest conservation.
#' @param target.temps Target melting temperatures for optimized primer sets in Celsius.
#' Only required when optimizing primers for both strand directions and one optimization was already performed.
#' @param required.cvg The target ratio of covered template sequences. 
#' If the target ratio cannot be reached, the constraint settings are relaxed up to the relaxation limits.
#' @param fw.primers List with optimized primer data frames corresponding to \code{target.temps}. 
#' Only required for optimizing both strand directions and only 
#' in the second optimization run in order to check for cross dimerization.
#' @param cur.results.loc Directory for storing results of the primer design procedure.
#' @param primer.df A data frame of evaluated primer candidates that can be optimized directly.
#' @param updateProgress Shiny progress callback function.
#' @return A list containing the results of the primer design procedure:
#' \describe{
#' \item{\code{opti}:}{A \code{Primers} object representing the set of optimized primers.}
#' \item{\code{all_results}:}{A list containing the optimal results for each sampled melting temperature range
#'                            in terms of a \code{Primers} object 
#'                            in case that the \code{melting_temp_diff} constraint was active. Otherwise,
#'                            \code{all_results} only has a single entry representing a primer set relating to an undefined melting temperature.}
#' \item{\code{used_constraints}:}{A \code{DesignSettings} object with the (adjusted) analysis settings.}
#' \item{\code{filtered}:}{A \code{Primers} object containing the primer candidates that passed the filtering procedure
#'                         and which gave rise to the final optimal set.}
#' }
#' @keywords internal
design_primers.single <- function(template.df, sample.name, 
    mode.directionality = c("fw", "rev"), 
    settings, timeout, opti.algo, allowed.region.definition, init.algo, 
    max.degen, conservation, target.temps, required.cvg, 
    fw.primers = NULL, 
    cur.results.loc = NULL, primer.df = NULL, updateProgress = NULL) {
   
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
   
    ####################
    # 'Unlist' settings
    #######################
    constraint.settings <- filters(settings)
    constraint.limits <- filterLimits(settings)
    PCR.settings <- PCR(settings)
    primer_conc <- PCR.settings$primer_concentration
    template_conc <- PCR.settings$template_concentration
    na_salt_conc <- PCR.settings$Na_concentration
    mg_salt_conc <- PCR.settings$Mg_concentration
    k_salt_conc <- PCR.settings$K_concentration
    tris_salt_conc <- PCR.settings$Tris_concentration
    other.settings <- conOptions(settings)
    allowed.mismatches <- other.settings$allowed_mismatches
    allowed.stop.codons <- other.settings$allowed_stop_codons
    allowed.other.binding.ratio <- other.settings$allowed_other_binding_ratio
    disallowed.mismatch.pos <- other.settings$disallowed_mismatch_pos
    ######################
    primer.lengths <- constraint.settings$primer_length
    opti.results.loc <- NULL
    filtering.results.loc <- NULL
    if (length(cur.results.loc) != 0) {
        # create some diagnostic folders
        cur.results.loc <- file.path(cur.results.loc, mode.directionality)
        filtering.results.loc <- file.path(cur.results.loc, "filtering")
        opti.results.loc <- file.path(cur.results.loc, "optimization", opti.algo)
        dir.create(cur.results.loc, showWarnings = FALSE, recursive = TRUE)
        dir.create(filtering.results.loc, showWarnings = FALSE, recursive = TRUE)
        dir.create(opti.results.loc, showWarnings = FALSE, recursive = TRUE)
    }
    # load primers/initialize primers
    used.filtering.constraints <- NULL
    annealing.temp <- NULL  # no pre-determined annealing temperature for filtering should be required
    filtered.result <- NULL
    used.settings <- settings
    if (length(primer.df) == 0) {
        # initialize and filter primer set
		message("# Phase 1: Initialization of primers")
        if (is.function(updateProgress)) {
            detail <- paste("Initialization_", mode.directionality, sep = "")
            updateProgress(1/3, detail, "inc")
        }
        primer.df <- initialize.primer.set(template.df, sample.name, primer.lengths, allowed.region.definition, 
            mode.directionality, init.algo, max.degen, conservation, cur.results.loc)
        # filter primers
		message("# Phase 2: Selection of primers")
        if (is.function(updateProgress)) {
            detail <- paste("Filtering_", mode.directionality, sep = "")
            updateProgress(1/3, detail, "inc")
        }
        filtered.result <- filter.primer.set.opti(primer.df, sample.name, 
            template.df, settings, mode.directionality, required.cvg, 
            filtering.results.loc, target.temps)
        filtered.df <- filtered.result$data
        used.settings <- filtered.result$used_settings # update the setings object
    } else {
        filtered.df <- primer.df
    }
    if (is.function(updateProgress)) {
        detail <- paste("Optimizing_", mode.directionality, sep = "")
        updateProgress(1/3, detail, "inc")
    }
    # ensure that melting temperature is available if melting_temp_diff constraint is active:
    if ("melting_temp_diff" %in% names(opti(used.settings)) && !"melting_temp" %in% colnames(filtered.df) && 
        Sys.which("melting-batch") != "") {
        # compute melting temperatures if they don't exist
        melting.temps <- compute.melting.temps(filtered.df, primer_conc, na_salt_conc, 
            mg_salt_conc, k_salt_conc, tris_salt_conc, mode.directionality)
        filtered.df <- cbind(filtered.df, melting.temps)
    }
	message("# Phase 3: Optimization of primers")
    if (opti.algo == "ILP") {
        # ILP OPTIMIZATION optimize with ILPs
        optimal.primers.data <- optimize.ILP(filtered.df, template.df, used.settings,
            primer_conc, template_conc, na_salt_conc, mg_salt_conc, 
            k_salt_conc, tris_salt_conc, allowed.mismatches, 
            allowed.other.binding.ratio, allowed.stop.codons, allowed.region.definition, 
            disallowed.mismatch.pos, target.temps, required.cvg, fw.primers = fw.primers, 
            diagnostic.location = opti.results.loc, timeout = timeout)
    } else {
        # GREEDY OPTIMIZATION
        optimal.primers.data <- select.primers.by.cvg(filtered.df, used.settings, 
            template.df, mode.directionality, required.cvg, constraint.settings$primer_coverage["primer_mismatches.max"], 
            primer_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, 
            template_conc, allowed.other.binding.ratio, 
            allowed.stop.codons, allowed.region.definition, 
            disallowed.mismatch.pos, target.temps, fw.primers)
    }
    opti.set <- optimal.primers.data$opti
    cvg.s <- paste0(round(get_cvg_ratio(opti.set, template.df) * 100, 2), "%")
    message("\to Selected ", nrow(opti.set), " primers (", cvg.s, " coverage).")
    # update melting_temp_diff & cross-dimer entries in all optimized sets
    for (i in seq_along(optimal.primers.data$all_results)) {
        cur.set <- optimal.primers.data$all_results[[i]]
        cur.set <- update.opti.results(cur.set, used.settings, template.df)
        optimal.primers.data$all_results[[i]] <- cur.set
    }
    opti.set <- update.opti.results(opti.set, used.settings, template.df)
    optimal.primers.data$opti <- opti.set
    # add filtering results to output
    optimal.primers.data$filtered <- filtered.result
    ###########
    # create a list of settings for the output
    used.settings <- list(optimal.primers.data$used_constraints)
    names(used.settings) <- mode.directionality
    optimal.primers.data$used_constraints <- used.settings
    if (length(opti.results.loc) != 0) {
        message("Storing results")
        write.out.primer.info(opti.results.loc, optimal.primers.data, mode.directionality, 
        used.settings, sample.name, template.df, max.degen)
        visualize.all.results(sample.name, filtering.results.loc, opti.results.loc, primer_conc, 
        na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, 
        settings, mode.directionality, used.settings, required.cvg)
    }
    return(optimal.primers.data)
}

#' Visualization of Design Results.
#'
#' Visualizes all results from designing primers.
#'
#' @param sample Identifier of the design run.
#' @param filtering.results.loc Location of filtering results.
#' @param opti.results.loc Location of optimization results.
#' @param primer_conc Primer concentration.
#' @param template_conc Template concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param settings The \code{DesignSettings} object.
#' @param mode.directionality Strand direction for which primers were designed.
#' @param used.settings A list with the used settings for optimization (fields "fw" and "rev").
#' @param required.cvg The required coverage.
#' @return Writes visualizations to files in.
#' @keywords internal
visualize.all.results <- function(sample, filtering.results.loc, opti.results.loc, 
    primer_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, 
    settings, mode.directionality, used.settings, required.cvg) {
    if (length(filtering.results.loc) == 0 || length(opti.results.loc) == 0) {
        return(NULL)  # not suppoed to plot anything
    }
    #### Plot filtering stats
    excluded.df <- read_primers(file.path(filtering.results.loc, paste(sample, "_excluded_primers.csv", 
        sep = "")))
    template.df <- read_templates(file.path(filtering.results.loc, paste(sample, "_filtered_templates.csv", 
        sep = "")))
    filtered.df <- read_primers(file.path(filtering.results.loc, paste(sample, 
        "_filtered_primers.csv", sep = "")))
    filtered.stats <- read.csv(file.path(filtering.results.loc, paste(sample, "_filtered_primers_stats.csv", 
        sep = "")), stringsAsFactors = FALSE)
    relax.stats <- file.path(filtering.results.loc, paste(sample, "_filtered_primers_stats_relax.csv", 
        sep = ""))
    if (file.size(relax.stats) > 3) {
        stats.relax <- read.csv(file.path(filtering.results.loc, paste(sample, "_filtered_primers_stats_relax.csv", 
            sep = "")), stringsAsFactors = FALSE)
    } else {
        stats.relax <- NULL
    }
    visualize.filtering.results(sample, filtering.results.loc, mode.directionality, 
        excluded.df, template.df, filtered.df, filtered.stats, stats.relax, primer_conc, 
        na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, 
        settings)
    # optimization plots: differentiate between different optimization methods
    optimal.primers <- try(read_primers(file.path(opti.results.loc, paste(sample, 
        "_FINAL_optimized_primers.csv", sep = ""))))
    if (class(optimal.primers) == "try-error") {
        return()  # no available results here
    }
    optimal.primers[is.na(optimal.primers)] <- ""
    template.df <- update_template_cvg(template.df, optimal.primers)  # update lex cvg 
    plot.all.cvg.info(sample, opti.results.loc, optimal.primers, template.df, mode.directionality, 
        "optimized", primer_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, 
        settings, required.cvg, used.settings) 
}
#' Visualization of Filtering Results.
#'
#' Visualizes the filtering results.
#'
#' @param sample Primer design run identifier.
#' @param results.loc Location where the filtering results are stored.
#' @param mode.directionality Design direction.
#' @param excluded.df Data frame with excluded primers.
#' @param template.df Template data frame.
#' @param filtered.df Primer data frame containing the primers that passed the constraints.
#' @param stats.relax Filtering statistics after relaxation.
#' @param primer_conc Primer concentration.
#' @param template_conc Template concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param settings A \code{DesignSettings} object.
#' @param required.cvg The required coverage.
#' @return Write-out of filtering results.
#' @keywords internal
visualize.filtering.results <- function(sample, results.loc, mode.directionality, 
    excluded.df, template.df, filtered.df, filtered.stats, stats.relax, primer_conc, na_salt_conc, 
    mg_salt_conc, k_salt_conc, tris_salt_conc, settings, required.cvg) {
    plot.all.filtering.stats(results.loc, sample, excluded.df, filtered.stats, stats.relax, 
        template.df)
    # create coverage plots/stats
    plot.all.cvg.info(sample, results.loc, filtered.df, template.df, mode.directionality, 
        "filtering", primer_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, 
        settings, required.cvg)
    cvg.stats <- get_cvg_stats(filtered.df, template.df)
    write.csv(cvg.stats, file.path(results.loc, paste(sample, "_filtering_cvg_stats.csv", 
        sep = "")), row.names = FALSE)
}

#' Plot Filtering Stats.
#'
#' Plots filtering statistics.
#'
#' @param results.loc Location where the filtering results are stored.
#' @param sample Primer design run identifier.
#' @param excluded.df Data frame with excluded primers.
#' @param filtered.stats Filtering statistics data frame.
#' @param stats.relax Filtering statistics after relaxation.
#' @param template.df Template data frame.
#' @return Write-out of filtering results.
#' @keywords internal
plot.all.filtering.stats <- function(results.loc, sample, excluded.df, filtered.stats, 
    stats.relax, template.df) {
    # Plot all filtering stats
    p <- plot.excluded.hist(excluded.df, filtered.stats, template.df)
    my_ggsave(file.path(results.loc, paste(sample, "_excluded_hist.png", sep = "")), p, 
        units = "cm", width = 40, height = 20)
    p <- plot.filtering.stats(filtered.stats, stats.relax)
    my_ggsave(file.path(results.loc, paste(sample, "_filtering_stats.png", sep = "")), p,
        units = "cm", width = 30, height = 15)
    p <- plot.filtering.stats.cvg(filtered.stats, stats.relax)
    my_ggsave(file.path(results.loc, paste(sample, "_filtering_stats_cvg.png", sep = "")), p,
        units = "cm", width = 30, height = 15)
    p <- plot.filtering.runtime(filtered.stats)
    my_ggsave(file.path(results.loc, paste(sample, "_filtering_runtime.png", sep = "")), p)
} 
#' Plots Coverage Information
#'
#' Visualizes all coverage-related data.
#'
#' 
#' @param sample Primer design run identifier.
#' @param results.loc Location where the filtering results are stored.
#' @param primers Primer data frame.
#' @param template.df Template data frame.
#' @param mode.directionality Design direction.
#' @param identifier Identifies whether filtering or optimization info should be displayed.
#' @param primer_conc Primer concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param settings The \code{DesignSettings} object.
#' @param required.cvg The required coverage.
#' @param used.settings A list containing \code{DesignSettings} objects for the 'fw' and 'rev' optimization.
#' @return Writes plots to files.
#' @keywords internal
plot.all.cvg.info <- function(sample, results.loc, primers, template.df, mode.directionality, 
    identifier = c("filtering", "optimized"), primer_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, 
    settings, required.cvg, used.settings = NULL) {

    if (length(identifier) == 0) {
        stop("Please supply the 'identifier' arg.")
    }
    identifier <- match.arg(identifier)
    if (length(primers) == 0 || length(template.df) == 0 || length(results.loc) == 0) {
        return(NULL)
    }
    constraint.settings <- constraints(settings)
    opti.constraints <- opti(settings)
    out.loc <- file.path(results.loc, "plots")
    dir.create(out.loc, showWarnings = FALSE)
    # primer group coverage overview table
    cvg.stats <- get_cvg_stats(primers, template.df)
    write.csv(cvg.stats, file = file.path(out.loc, paste(sample, "_", identifier, 
        "_cvg_stats.csv", sep = "")), row.names = FALSE)
    plot_template_cvg(primers, template.df)
    my_ggsave(file.path(out.loc, paste(sample, "_", identifier, "_template_cvg.png", sep = "")), 
        units = "cm", width = 15, height = 15)
    plot_template_cvg(primers, template.df, per.mismatch = TRUE)
    my_ggsave(file.path(out.loc, paste(sample, "_", identifier, "_template_cvg_mismatch.png", sep = "")), 
        units = "cm", width = 15, height = 15)
    # primer binding regions
    if (mode.directionality == "both") {
        relations <- c("fw", "rev")  # forward primer info
    } else {
        relations <- mode.directionality
    }
    for (i in seq_along(relations)) {
        relation <- relations[i]
        plot_primer_binding_regions(primers, template.df, relation, group = "all")
        my_ggsave(file.path(out.loc, paste(sample, "_", identifier, "_primer_binding_regions_", relation, ".png", 
            sep = "")), units = "cm", width = 30, height = 15)
    }
    if (identifier != "filtering") {
        # don't plot for the filtered primers: too many ..
        # create pdf report
        file <- file.path(out.loc, get_report_fname("report", sample))
        create_report(primers, template.df, file, settings, 
                    sample.name = sample, used.settings = used.settings, 
                    required.cvg = required.cvg)
        plot_primer_cvg(primers, template.df)
        my_ggsave(file.path(out.loc, paste(sample, "_", identifier, "_primer_cvg.png", 
            sep = "")), units = "cm", width = 25, height = 15)
        for (i in seq_along(relations)) {
            relation <- relations[i]
            p <- plot_primer(primers, template.df, primers$Identifier, relation)
            my_ggsave(file.path(out.loc, paste(sample, "_", identifier, "_", relation, "_primer_plot.png", 
                sep = "")), plot = p, height = 90, units = "cm", width = 30, dpi = 300)
        }
        # self dimerization
        self.dimer.cutoff <- ifelse("self_dimerization" %in% names(constraint.settings), 
            constraint.settings$self_dimerization, NA)
        annealing.temp <- compute_annealing_temp(primers, mode.directionality, 
                                template.df, PCR(settings)$Na_concentration,
                                PCR(settings)$Mg_concentration, 
                                PCR(settings)$K_concentration, 
                                PCR(settings)$Tris_concentration, 
                                PCR(settings)$primer_concentration)
        if (!is.na(self.dimer.cutoff)) {

            # only plot if data were supposed to be computed
            self.dimer.data <- compute.all.self.dimers(primers, primer_conc, 
                na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, annealing.temp)
            se <- self.dimer.data
            #se <- view.self.dimer.table(self.dimer.data)
            write.csv(se, file = file.path(out.loc, paste(sample, "_", identifier, "_self_dimer_details.csv", 
                sep = "")), row.names = FALSE)
            write.csv(self.dimer.data, file = file.path(out.loc, paste(sample, "_", identifier, 
                "_self_dimer_details_raw.csv", sep = "")), row.names = FALSE)
            p <- plot.dimer.dist(se, self.dimer.cutoff)
            my_ggsave(file.path(out.loc, paste(sample, "_", identifier, "_self_dimer_DeltaG.png", 
                sep = "")), plot = p)
            self.dim.table <- dimerization.table(se, self.dimer.cutoff, "Self-Dimerization")
            write.csv(self.dim.table, file = file.path(out.loc, paste(sample, "_", identifier, 
                "_self_dimer_counts.csv", sep = "")), row.names = FALSE)
            # cross dimerization
            cross.dimer.cutoff <- ifelse("cross_dimerization" %in% names(opti.constraints),
                opti.constraints$cross_dimerization, -8)
            all.cross.dimer.data <- compute.all.cross.dimers.unfiltered(primers, primer_conc, 
                na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, annealing.temp)
            cross.dimer.data <- compute.all.cross.dimers(primers, primer_conc, 
                na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, annealing.temp,
                all.cross.dimer.data)
            if (length(cross.dimer.data) != 0) {
                cross.dimer.data <- na.omit(cross.dimer.data)
            }
            se <- cross.dimer.data
            #se <- view.cross.dimer.table(cross.dimer.data)
            write.csv(se, file = file.path(out.loc, paste(sample, "_", identifier, "_cross_dimer_details.csv", 
                sep = "")), row.names = FALSE)
            write.csv(cross.dimer.data, file = file.path(out.loc, paste(sample, "_", 
                identifier, "_cross_dimer_details_raw.csv", sep = "")), row.names = FALSE)
            p <- plot.dimer.dist(se, cross.dimer.cutoff)
            my_ggsave(file.path(out.loc, paste(sample, "_", identifier, "_cross_dimer_DeltaG.png", 
                sep = "")), plot = p)
            cross.dim.table <- dimerization.table(se, cross.dimer.cutoff, "Cross-Dimerization")
            write.csv(cross.dim.table, file = file.path(out.loc, paste(sample, "_", identifier, 
                "_cross_dimer_counts.csv", sep = "")), row.names = FALSE)
        }
        # mismatches
        mm.table.fw <- compute.mismatch.table(primers, template.df, "fw")
        write.csv(mm.table.fw, file = file.path(out.loc, paste(sample, "_", identifier, 
            "_mismatch_priming.csv", sep = "")), row.names = FALSE)
        # constraint evaluation: which constraints are fulfilled?
        p <- plot_constraint_fulfillment(primers, settings)
        my_ggsave(file.path(out.loc, paste(sample, "_", identifier, "_constraints_evaluation.png", 
            sep = "")), plot = p)
        # output optimal set subsets
        groups <- NULL  # consider cvg of all groups
        compute.all.primer.subsets.ILP(primers, template.df, 1, groups, out.loc, required.cvg)
    }
}
#' Computation of Primer Subsets
#'
#' Computes all optimal primer subsets and stores their plots.
#'
#' @param primer.df Primer data frame.
#' @param template.df Template data frame.
#' @param k Subset size-increment.
#' @param groups Identifiers of template groups in order to limit coverage to
#' certain groups of template sequences.
#' @param cur.results.loc Location for storing the results.
#' @param required.cvg The required coverage ratio.
#' @return Write-out of results.
#' @keywords internal
compute.all.primer.subsets.ILP <- function(primer.df, template.df, k, groups, cur.results.loc, required.cvg = 1) {
    if (length(primer.df) == 0 || nrow(primer.df) == 0) {
        return(NULL)
    }
    identifier <- "subsets"
    if (length(groups) != 0) {
        identifier <- paste(identifier, "_groups_", paste(groups, collapse = "_", 
            sep = ""), sep = "")
    }
    top.primers.overall <- subset_primer_set(primer.df, template.df, k, groups, 
                                            identifier, cur.results.loc)
    p <- plot_primer_subsets(top.primers.overall, template.df, 
                        required.cvg = required.cvg)
    my_ggsave(file.path(cur.results.loc, paste(identifier, "_coverage_ramp.png", sep = "")), p,
        units = "cm", width = 15, height = 15)
}
#' Creation of Melting Temperature Groups
#'
#' Creates a data frame identifying target melting temperatures of individual primer sets.
#'
#' @param primers An object of class \code{Primers} for which to create groups based on melting temperatures.
#' @param template.df An object of class \code{Templates} corresponding to the primers.
#' @param settings A \code{DesignSettings} objects.
#' @param target.temps Pre-defined target melting temperatures to use instead
#' of automatically determining groups from the \code{primers}.
#' @return Data frame with target melting temperatures for individual primer sets.
#' @keywords internal
create.Tm.brackets <- function(primers, template.df, settings, target.temps = NULL) {
    group.assignment <- rep(1, nrow(primers)) # all in one group by default
    group.df <- NULL
    group.df <- data.frame(GroupID = 1, min_Tm = NA, max_Tm = NA, 
                            target_Tm = NA, annealing_temp = NA)
    opti.constraints <- opti(settings)
    max.Tm.delta <- opti.constraints[["melting_temp_diff"]]["max"]
    # ensure that melting temperatures are available if melting_temp_diff is to be constrained
    if ("melting_temp_diff" %in% names(opti.constraints) && !"melting_temp" %in% colnames(primers) && nrow(primers) != 0) {
        # melting temp is not available although melting_temp_diff should be respected -> compute melting temperatures here
        primers <- check_constraints(primers, template.df, settings, active.constraints = "melting_temp_diff")
    } else if ("melting_temp" %in% colnames(primers) && "melting_temp_diff" %in% names(opti.constraints)) {
        # we have a melting temp, but we might still have missing entries ...
        na.idx <- which(is.na(primers$melting_temp))
        if (length(na.idx) != 0) {
            primers[na.idx,] <- check_constraints(primers[na.idx,], template.df, settings, active.constraints = "melting_temp_diff")
        }
    }
    if (length(target.temps) != 0 && "melting_temp_diff" %in% names(opti.constraints)) {
        # init group.df with target temps
        group.df <- data.frame(GroupID = seq_along(target.temps), min_Tm = target.temps - (max.Tm.delta/2), 
                                max_Tm = target.temps + (max.Tm.delta/2), target_Tm = target.temps, annealing_temp = NA)
    } else if ("melting_temp_diff" %in% names(opti.constraints) && nrow(primers) != 0) {
        # no target temperatures specified -> determine grouping automatically:
        # ensure that groups overlap by selecting steps half the size of the allowed temperature difference
        groups <- c(seq(min(primers$melting_temp, na.rm = TRUE), max(primers$melting_temp, na.rm = TRUE), ifelse(max.Tm.delta !=0, max.Tm.delta / 2, 1)), max(primers$melting_temp, na.rm = TRUE))
        if (length(groups) == 2) {
            # ensure that 'group.df' can be constructed
            groups <- c(groups, max(primers$melting_temp, na.rm = TRUE))
        }
        group.df <- data.frame(GroupID = 1:(length(groups) - 2), min_Tm = groups[1:(length(groups) - 
            2)], max_Tm = groups[3:length(groups)])
        group.df$target_Tm <- group.df$min_Tm + (group.df$max_Tm - group.df$min_Tm)/2
    }
    # assign annealing temperatures
    if (nrow(primers) != 0 && "melting_temp_diff" %in% names(opti.constraints)) {
        group.df$annealing_temp <- sapply(group.df$min_Tm, function(x) annealing.temp.rule.of.thumb(x))
        # assign primers to groups if melting temperature is available
        group.assignment <- sapply(primers$melting_temp, function(x) group.df$GroupID[which(x >= 
            group.df$min_Tm & x < group.df$max_Tm)])
        group.assignment <- sapply(group.assignment, function(x) ifelse(length(x) == 
            0, NA, x))
        mean.temp <- mean(primers$melting_temp)  # mean melting temperature of the data set serves as the starting position of temperature optimization
        target.temp.order <- order(abs(group.df$target_Tm - mean.temp))
        group.df <- group.df[target.temp.order, ]
    } else if (nrow(primers) != 0 && !"melting_temp_diff"  %in% names(opti.constraints)) {
        # all primers are in one group -> compute annealing temp for all
        # necessary for computations that depend on the Ta
        mode.directionality <- get.analysis.mode(primers)
        annealing.temp <- compute_annealing_temp(primers, mode.directionality, 
                                template.df, PCR(settings)$Na_concentration,
                                PCR(settings)$Mg_concentration, 
                                PCR(settings)$K_concentration, 
                                PCR(settings)$Tris_concentration, 
                                PCR(settings)$primer_concentration)
        group.df[, "annealing_temp"] <- min(annealing.temp)
        # the melting temp range doesn't need to be computed in this case
        # no adjustment of groups necessary
    }
    #print("Tm brackets: primer Tms:")
    #print(primers$melting_temp) # check whether all have a melting temperature
    return(list(group_idx = group.assignment, df = group.df, primers = primers))
}
#' Cross-Dimerization Filtering
#'
#' Removes cross-dimerizing primers from the input data.
#'
#' @param primers.rev The primer data set to be filtered for cross-dimers.
#' @param primers.fw The primer data set to be checked against cross-dimerization.
#' @param opti.constraints List with optimization constraint settings.
#' @param annealing.temp The PCR annealing temperature.
#' @param primer_conc Primer concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @return \code{primers.rev} with removed cross-dimerizing primers.
#' @keywords internal
#filter.cross.comp <- function(primers.rev, primers.fw, opti.constraints, annealing.temp, 
    #primer_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc) {
    ## determines whether primers in first set (rev) dimerize with primers in second
    ## set (fw) and removes primers from first set. this function is used during
    ## optimization, when one direction has already been optimized (primers.fw) and
    ## the other one hasn't been yet.  no adjustment necessary ...
    #if (!"cross_dimerization" %in% names(opti.constraints) || length(primers.fw) == 
        #0 || length(primers.rev) == 0) {
        #return(primers.rev)
    #}
    #primers.1 <- primers.rev$Reverse  # rows: selection
    #primers.2 <- primers.fw$Forward  # cols
    #G.df <- get.cross.dimers(primers.1, primers.2, primer_conc, na_salt_conc, 
        #mg_salt_conc, k_salt_conc, tris_salt_conc, annealing.temp, mode = "asymmetric")
    #deltaG.cutoff <- opti.constraints$cross_dimerization["min"]
    ## problem with create.G.matrix: dimension is defined by primer.df ... always
    # square!
    #G.matrix <- create.G.matrix(primers.rev, G.df, primers.fw)  # min deltaG values of cross-dimerization conformations for every pair of primers
    #D <- compute.dimer.matrix(G.matrix, deltaG.cutoff)  # 1 for dimers, 0 else
    ## filter:
    #rm.idx <- which(apply(D, 1, function(x) any(x == 1)))
    #if (length(rm.idx) != 0) {
        #primers.rev <- primers.rev[-rm.idx, ]
    #}
    #return(primers.rev)
#}
#' Creation of Melting Temperature Sets.
#'
#' Stratifies primers according to their melting temperatures and checks temperature-dependent constraints.
#'
#' Cross-dimerization between forward and reverse primers is checked here, in case this is a second 
#' optimization run and a list \code{primers.fw} with primer data frames is provided.
#'
#' @param primer.df Primer data frame.
#' @param template.df Template data frame.
#' @param Tm.brackets Data frame with target primer melting temperatures.
#' @param settings A \code{DesignSettings} object.
#' @param mode.directionality Identifier of strand for which primers shall be designed.
#' @param primer_conc Primer concentration.
#' @param template_conc Template concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param allowed.mismatches The number of mismatches primers are allowed to have with the templates.
#' @param allowed.other.binding.ratio Ratio of primers allowed to bind to non-target regions.
#' @param allowed.stop.codons  Consider mismatch binding events that induce stop codons.
#' @param allowed.region.definition Definition of the allowed region.
#' @param disallowed.mismatch.pos The number of positions from the primer 3' end where mismatches should not be allowed.
#' All primers binding templates with mismatches within \code{disallowed.mismatch.pos} from the 3' end are disregarded.
#' @param opti.mode Compute optimization constraints and relax delta Tm if necessary.
#' @param required.cvg Target coverage ratio.
#' @param primers.fw Already designed primer sets for the target temperatures given in \code{Tm.brackets}.
#' Used to determine cross-dimerization.
#' @param diagnostic.location Directory for storing results.
#' @param updateProgress Shiny progress callback function.
#' @return Primer data frames for every target temperature.
#' @keywords internal
compute.Tm.sets <- function(primer.df, template.df, Tm.brackets, settings, mode.directionality = c("fw", "rev"), 
    primer_conc, template_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, 
    allowed.mismatches, allowed.other.binding.ratio, 
    allowed.stop.codons, allowed.region.definition, disallowed.mismatch.pos, 
    opti.mode = FALSE, required.cvg = NULL, 
    primers.fw = NULL, diagnostic.location = NULL, updateProgress = NULL) {
  
    # determine constraints that need to be computed: consider only primer efficiency
    # and secondary struct from opti constraints
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' arg.")
    }
    mode.directionality <- match.arg(mode.directionality) 
    # load initial settings:
    opti.constraints <- opti(settings)
    opti.limits <- optiLimits(settings)
    comp.constraints <- NULL  # active constraints for filtering during optimization (annealing-temp dependent)
    if (opti.mode) {
        # compute constraint values for selected constraints
        if ("secondary_structure" %in% names(opti.constraints)) {
            # should be computed in filters now ...
            comp.constraints <- c(comp.constraints, "secondary_structure")
        }
    }
    target.temps <- Tm.brackets$df$target_Tm
    Tm.sets <- vector("list", length(target.temps))
    out.settings <- vector("list", length(target.temps))
    if (opti.mode) {
        message("Determining melting temperature sets:")
    } else {
        message("Checking melting temperature set coverage:")
    }
    for (i in seq_along(target.temps)) {
        # don't parallelize here, it's more important to parallelize constraint computations
        cur.opti.constraints <- opti.constraints
        cur.opti.limits <- opti.limits
        #cat(paste('Opti constraint filtering\n\to Iteration: ', i, '/',
            #length(target.temps), '\n\to Target temperature: ', target.temps[i], 
            #\n\to Constraints: ', paste(comp.constraints, collapse = ','), '\n', sep= ''))
        target.temp <- target.temps[i]
        # compute absolute melting temperature differences
        annealing.temp <- Tm.brackets$df$annealing_temp[i]
        # select Tm set:
        if (!is.na(target.temp) & "melting_temp" %in% colnames(primer.df)) {
            # target is not specified or melting temperature is not available
            cur.sel <- which(abs(primer.df$melting_temp - target.temp) <= (cur.opti.constraints$melting_temp_diff["max"]/2))  # was based on opti limits before, but this was removed so that ILP does not need temperature constraint
        } else {
            cur.sel <- seq_len(nrow(primer.df))  # select all
        }
         # message(paste('Available primers in Tm set:', length(cur.sel), sep = '')) no
        # available primers for current temp
        Tm.set <- primer.df[cur.sel, ]
        message("Melting temperature set @ ", target.temps[i], " with ",
                 nrow(Tm.set), " primers (", get_cvg_ratio(Tm.set, template.df, as.char = TRUE), " coverage)")
        if (length(Tm.set) == 0 || nrow(Tm.set) == 0) {
            # set to empty data frame!
            Tm.sets[[i]] <- Tm.set
            out.settings[[i]] <- settings
            next
        }
        if (opti.mode) {
            ###########
            # evaluate and relax constraints
            ############
            # try not to lose any coverage:
            # since we may have just selected 1 set with 'required.cvg' we cannot
            # try to maintain the existing cvg -> at least one set provided by the relaxation should achieve the specified target cvg from the design function
            target.cvg <- min(get_cvg_ratio(primer.df, template.df), required.cvg) # try to attain the required cvg if possible, or less, if we cannot anymore.
            message("\to Computing optimization constraints for target coverage ", paste0(round(target.cvg * 100, 2), "%"))
            # 1. Melting temp difference relaxation:
            while (get_cvg_ratio(Tm.set, template.df) < target.cvg && "melting_temp_diff" %in% names(opti.constraints) && "melting_temp" %in% colnames(primer.df)) {
                # relax Delta Tm
                cur.opti.constraints$melting_temp_diff <- relax.opti.constraints(cur.opti.constraints, opti.limits, opti.constraints)$melting_temp_diff
                cur.opti.limits$melting_temp_diff <- set.new.limits(cur.opti.limits, opti.limits, opti.constraints)$melting_temp_diff

                # create new Tm set according to relaxed DeltaG setting:
                cur.sel <- which(abs(primer.df$melting_temp - target.temp) <= (cur.opti.constraints$melting_temp_diff["max"]/2))
                Tm.set <- primer.df[cur.sel,]
            }
            #message("After melting temp deviation filter: ", get_cvg_ratio(Tm.set, template.df))
            # 2. Computation of cross-terms between fw/rev primers
            if (mode.directionality == "rev" && length(primers.fw) != 0 && 
                    "cross_dimerization" %in% names(cur.opti.constraints)) {
                # consider cross-terms of fw and rev primers when optimizing both fw and rev primers
                ions <- compute.sodium.equivalent.conc(na_salt_conc, mg_salt_conc, 
                                           k_salt_conc, tris_salt_conc)
                G.df <- get.cross.dimers(Tm.set$Reverse, primers.fw[[i]]$Forward, 
                                ions, annealing.temp, 
                                no.structures = TRUE, mode = "asymmetric")
                G.matrix <- create.G.matrix(Tm.set, G.df, primers.fw[[i]])  # min deltaG values of cross-dimerization conformations for every pair of primers
                D <- compute.dimer.matrix(G.matrix, cur.opti.constraints$cross_dimerization)  # 1 for dimers, 0 else
                # filr new Tm.set -> check cvg
                rm.idx <- which(apply(D, 1, function(x) any(x == 1)))
                new.Tm.set <- Tm.set
                if (length(rm.idx) != 0) {
                    new.Tm.set <- new.Tm.set[-rm.idx, ]
                }
                # relax cross term conditions:
                while (get_cvg_ratio(new.Tm.set, template.df) < target.cvg) { # relax and re-evaluate
                    cur.opti.constraints$cross_dimerization <- relax.opti.constraints(cur.opti.constraints, opti.limits, opti.constraints)$cross_dimerization
                    cur.opti.limits$cross_dimerization <- set.new.limits(cur.opti.limits, opti.limits, opti.constraints)$cross_dimerization
                    D <- compute.dimer.matrix(G.matrix, cur.opti.constraints$cross_dimerization)
                    rm.idx <- which(apply(D, 1, function(x) any(x == 1)))
                    new.Tm.set <- Tm.set
                    if (length(rm.idx) != 0) {
                        new.Tm.set <- new.Tm.set[-rm.idx, ]
                    }
                }
                Tm.set <- new.Tm.set # set to filtered df
            }
            #message("After cross dimerization filter: ", get_cvg_ratio(Tm.set, template.df))
            # 3. Computation of other Constraints, e.g. primer efficiency
            # there's no relaxation here. the filtered data is only used if we can still reach the target coverage.
            constraint.df <- compute.constraints(Tm.set, mode.directionality, template.df, settings, updateProgress)
            Tm.set.data <- filter.by.constraints(Tm.set, constraint.df, cur.opti.constraints[comp.constraints], 
                comp.constraints, mode.directionality, template.df)
            new.Tm.set <- Tm.set.data$Filtered
            if (get_cvg_ratio(new.Tm.set, template.df) >= target.cvg ) {
                # only use efficiency-filtered data, if we're at the target cvg, otherwise use the non-filtered data
                Tm.set <- new.Tm.set
            }
            #message("After cross dimerization filter: ", get_cvg_ratio(Tm.set, template.df))
            # write out filtered data set to diagnostic location output info
            if (length(diagnostic.location) != 0) {
                write.csv(Tm.set, file = file.path(diagnostic.location, paste("filtered_opti_data_target_temp_", 
                    target.temp, ".csv", sep = "")), row.names = FALSE)  # write every iteration in case we have to stop the analysis due to time constraints
                write.csv(Tm.set.data$Excluded, file = file.path(diagnostic.location, 
                    paste("filtered_opti_data_excluded_target_temp_", target.temp, ".csv", 
                      sep = "")), row.names = FALSE)  # write every iteration in case we have to stop the analysis due to time constraints
            }
        }
        # output relaxed settings and sets
        cur.settings <- settings
        constraintLimits(cur.settings)[names(cur.opti.limits)] <- cur.opti.limits
        constraints(cur.settings)[names(cur.opti.constraints)] <- cur.opti.constraints
        out.settings[[i]] <- cur.settings
        Tm.sets[[i]] <- Tm.set # store melting temperature-target set
    }
    result <- list("sets" = Tm.sets, "settings" = out.settings)
    return(result)
}
#' Relaxation of Optimization Constraints
#'
#' Relax optimization constraints.
#'
#' @param cur.opti.constraints List with optimization constraint settings.
#' @param initial.opti.limits Initial optimization limits.
#' @param initial.opti.constraints Initial optimization constraints.
#' @return Relaxed optimization constraints.
#' @keywords internal
relax.opti.constraints <- function(cur.opti.constraints, initial.opti.limits, initial.opti.constraints) {
    if (length(cur.opti.constraints) == 0) {
        # nothing to relax ..
        return(NULL)
    }
    cur.opti.constraints <- set.new.limits(cur.opti.constraints, initial.opti.limits, initial.opti.constraints)
    return(cur.opti.constraints)
}

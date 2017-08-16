
penalize_primer <- function(deviations, alpha = 0.5) {
    # alpha: tuning parameter for penalizing maximum deviation
    if (!is.numeric(alpha) || alpha <0 || alpha > 1) {
        stop("Please provide a numeric in the interval [0,1] for 'alpha'.")
    }
    abs.dev <- abs(deviations)
    max.idx <- which.max(abs.dev)
    # minimal penalty is 0. 
    score <- (alpha * abs.dev[max.idx]) + (1 - alpha) * sum(abs.dev)
    return(score)
}
#' Scoring of Primers.
#'
#' Computes scores for a set of primers based on the deviations
#' of the primers from the constraints.
#'
#' The penalty of a primer is computed in the following way. Let 
#' \code{d} be a vector indicating the absolute deviations from individual constraints
#' and let \code{p} be the scalar penalty that is assigned to a primer. We define
#' \deqn{p = \alpha \cdot \max_i d_i + \sum_i (1 - \alpha) \cdot d_i}
#' such that for large values of \code{alpha} the maximal deviation 
#' dominates giving rise to a local penalty (reflecting the largest absolute deviation) and for small 
#' \code{alpha} the total deviation dominates giving rise to a global penalty 
#' (reflecting the sum of constraint deviations). 
#' When \code{alpha} is 1 only the most extreme absolute deviation is considered and
#' when \code{alpha} is 0 the sum of all absolute deviations is computed.
#' @param primer.df A \code{Primers} object containing the primers
#' that are to be scored.
#' @param settings A \code{DesignSettings} object containing the
#' settings that are evaluated when computing the deviation.
#' @param active.constraints A character vector of constraint identifiers
#' that are considered for scoring the primers.
#' @param alpha A numeric that is used to determine the trade-off
#' between the impact of the maximal observed deviation and the total
#' deviation. At its default \code{alpha} is set to 0.5 such that
#' the maximal deviation and the total deviation have an equal weight
#' when computing the penalties.
#' @return A data frame containing scores for the primers.
#' @export
#' @examples
#' data(Ippolito)
#' primer.scores <- score_primers(primer.df, settings)
score_primers <- function(primer.df, settings, active.constraints = names(constraints(settings)), alpha = 0.5) {
    if (!is(primer.df, "Primers")) {
        stop("Please input a 'Primers' object as 'primer.df'.")
    }
    my_penalty <- function(deviation) {
        # local function where 'alpha' doesn't need to be supplied for ddply
        penalize_primer(deviation, alpha = alpha)
    }
    constraint.settings <- constraints(settings)[active.constraints]
    deviations <- get_constraint_deviation_data(primer.df, constraint.settings)
    penalty.df <- plyr::ddply(deviations, c("ID"), 
                              plyr::here(plyr::summarize), Penalty = my_penalty(substitute(Deviation)), Deviation = sum(abs(substitute(Deviation))))
    return(penalty.df)
}
design_primers_global_single <- function(template.df, settings, mode.directionality = c("fw", "rev"), 
                                         sample.name = NULL, target.temps = NULL, init.algo = c("naive", "tree"),
                                         max.degen = 16, conservation = 1, cur.results.loc = NULL, prefilter = FALSE,
                                         max.penalties = c(0, 0.25, 0.5, 0.75, 1, 1.5, 2, 3, 5, Inf)) {
    # load example primers:
    if (length(sample.name) == 0) {
        sample.name <- "TEST"
    }
    cur.results.loc <- NULL
    primer.df <- initialize.primer.set(template.df, sample.name, constraints(settings)$primer_length,  # N = 10,707
                conOptions(settings)$allowed_region_definition, mode.directionality, init.algo, max.degen, conservation, cur.results.loc)
    message("Number of primers: ", nrow(primer.df))
    # don't sub-sample the primer set
    #my.sel <- sample(seq_len(nrow(primer.df)), 50)
    #primer.df <- primer.df[my.sel, ]
    # 0) Pre-filter? Necessary to reduce possible number of primers sometimes (> 10, 000 makes problems)
    if (prefilter) {
        filter.settings <- settings
        filter.constraints <- list("gc_ratio" = c("min" = 0.3, "max" = 0.7), # ~10%
                                "melting_temp_range" = c("min" = 40, "max" = 67), # ~3%
                                "no_runs" = c("min" = 0, "max" = 5), # ~0%
                                "no_repeats" = c("min" = 0, "max" = 5)) # ~0%
        constraints(filter.settings) <- filter.constraints
        primer.df <- filter_primers(primer.df, template.df, filter.settings)
        message("Number of primers after filtering: ", nrow(primer.df)) # N = 9,930 > ok? -> deal with cross dimer computation? .. could use greedy? or filter more
    }
    # A) Determine all constraints for all primers
    # select constraints that need to be computed:
    exclude.constraints.check <- c("primer_length", "gc_ratio", 
                             "no_runs", "no_repeats",
                             "melting_temp_diff",
                             "cross_dimerization")
    compute.constraints <- setdiff(names(constraints(settings)), exclude.constraints.check)
    primer.df <- check_constraints(primer.df, template.df, settings, 
                                   compute.constraints)
    message("Done with evaluating the primers!")
    # B) Score all primers
    exclude.constraints <- c("primer_length", "gc_ratio", 
                             "no_runs", "no_repeats",
                             "melting_temp_range",
                             "primer_coverage", "melting_temp_diff",
                             "cross_dimerization")
    score.constraints <- setdiff(names(constraints(settings)), exclude.constraints)
    penalty.df <- score_primers(primer.df, settings, active.constraints = score.constraints)
    penalties <- penalty.df$Penalty[match(primer.df$ID, penalty.df$ID)]
    primer.df$Penalty <- penalties
    if (length(cur.results.loc != 0)) {
        png(file.path(cur.results.loc, paste0("design_penalties_", sample.name, ".png")))
        hist(penalties)
        dev.off()
        save(primer.df, file = file.path(cur.results.loc, paste0("eval_primers_", sample.name, ".Rdata")))
    }
    #load(file.path("data", "global_primers.Rdata"))
    # C) Optimize (NB: consider melting_temp_diff, cross_dimerization)!
    # build dimerization matrix
    Tm.brackets <- create.Tm.brackets(primer.df, template.df, settings, target.temps)
    primer_conc <- PCR(settings)$primer_concentration
    na_salt_conc <- PCR(settings)$Na_concentration
    mg_salt_conc <- PCR(settings)$Mg_concentration
    k_salt_conc <- PCR(settings)$K_concentration
    tris_salt_conc <- PCR(settings)$Tris_concentration
    annealing.temp <- min(Tm.brackets$df$annealing_temp) # conservative dimerization calls
    message("Computing cross dimers ...")
    G.df <- compute.all.cross.dimers(primer.df, primer_conc, na_salt_conc, 
                            mg_salt_conc, k_salt_conc, tris_salt_conc, annealing.temp,
                            no.structures = TRUE)
    if (length(cur.results.loc) != 0) {
        save(G.df, file = file.path("data", paste0("dimer_energies_", sample.name, ".Rdata")))
    }
    #load(file.path("data", "global_dimer_energy.Rdata"))
    G.matrix <- create.G.matrix(primer.df, G.df) 
    message("Cross dimers computed!")
    cvg.matrix <- get.coverage.matrix(primer.df, template.df) # TODO: this could be a bit faster ..
    max.temp.diff <- constraints(settings)$melting_temp_diff["max"]
    dimer.cutoff <- constraints(settings)$cross_dimerization["min"]
    penalty.sets <- vector("list", length(max.penalties))
    message("Optimizing ...")
    for (i in seq_along(max.penalties)) {
        print(i)
        # explore several maximal penalties:
        max.penalty <- max.penalties[i]
        if (is.infinite(max.penalty)) {
            cur.max.temp.diff <- 50
        } else {
            cur.max.temp.diff <- max.temp.diff + (max.penalty * max.temp.diff)
        }
        cur.settings <- settings
        constraints(cur.settings)$melting_temp_diff <- c("max" = unname(cur.max.temp.diff))
        cur.temp.brackets <- create.Tm.brackets(primer.df, template.df, cur.settings, target.temps)
        cur.dimer.cutoff <- dimer.cutoff + (max.penalty * dimer.cutoff)
        D <- compute.dimer.matrix(G.matrix, cur.dimer.cutoff)
        temperature.sets <- vector("list", length(nrow(cur.temp.brackets))) 
        # define a temperature range
        for (j in seq_along(cur.temp.brackets$df$annealing_temp)) {
            # select only primers with penalties <= max.penalty and within melting temperature range
            target.Ta <- cur.temp.brackets$df$annealing_temp[j]
            run.id <- paste0("pen=", round(max.penalty, 2), "_Ta=", round(target.Ta,0))
            min.Tm <- cur.temp.brackets$df$min_Tm[j]
            max.Tm <- cur.temp.brackets$df$max_Tm[j]
            sel <- which(primer.df$Penalty <= max.penalty & primer.df$melting_temp <= max.Tm & primer.df$melting_temp >= min.Tm)
            if (length(sel) == 0) {
                # No primers available fulfilling the conditions
                next
            }
            cur.primers <- primer.df[sel,] 
            cur.D <- D[sel,sel, drop = FALSE]
            cur.cvg.matrix <- cvg.matrix[, sel, drop = FALSE]
            ILP <- ILPConstrained(cur.D, cur.cvg.matrix)
            #write.lp(ILP, "test_ILP.txt")
            return.val <- solve(ILP)
            if (return.val != 0) {
                # Optimization didn't work
                warning(return.val)
            }
            vars <- get.variables(ILP)
            primer.vars <- vars[seq_len(nrow(cur.primers))]
            sel.idx <- which(primer.vars == 1)
            opt.primers <- cur.primers[sel.idx, ]
            # set Run identifier
            if (length(opt.primers) != 0 && nrow(opt.primers) != 0) {
                opt.primers$Run <- run.id
            }
            print(as.numeric(get_cvg_ratio(opt.primers, template.df)))
            temperature.sets[[j]] <- opt.primers
        }
        penalty.sets[[i]] <- temperature.sets
    }
    message("Storing optimization results")
    primer.data <- unlist(penalty.sets, recursive = FALSE)
    # select non-empty sets only:
    sel <- which(unlist(lapply(primer.data, function(x) length(x) != 0 && nrow(x) != 0)))
    primer.data <- primer.data[sel]
    set.names <- unlist(lapply(primer.data, function(x) x$Run[1]))
    names(primer.data) <- set.names
    # enrich sets with missing constraint info:
    for (i in seq_along(primer.data)) {
        print(i)
        primer.data[[i]] <- check_constraints(primer.data[[i]], template.df, settings, exclude.constraints.check)
    }
    return(primer.data)
}

analyze_sets <- function() {
    mode.directionality <- "fw"
    target.temps <- NULL
    init.algo <- "tree"
    max.degen <- 16
    conservation <- 1
    sample.name <- "IGH"
    template.file <- system.file("extdata", "IMGT_data", "templates", "Homo_sapiens_IGH_functional_exon_exp.fasta", package = "openPrimeR")
    template.df <- read_templates(template.file, id.col = "GROUP", hdr.structure = c("ACCESSION", "GROUP"), delim = "|")
    # use the first 30 bases in the leader for binding
    analysis.type <- "leader"
    prefilter <- FALSE
    if (analysis.type == "first_30") {
        template.df <- assign_binding_regions(template.df, fw = c(1,30))
    } else if(analysis.type == "leader") {
        leader.loc <- system.file("extdata", "IMGT_data", "templates", "Homo_sapiens_IGH_functional_leader_exp.fasta", package = "openPrimeR")
        template.df <- assign_binding_regions(template.df, fw = leader.loc, rev = NULL)
        prefilter <- TRUE
    }
    # adjust settings:
    setting.file <- system.file("extdata", "settings", 
            "A_Taq_PCR_design.xml", package = "openPrimeR")
    settings <- read_settings(setting.file)
    # prevent amino acid mutations and stop codons
    cvg_constraints(settings)$substitution["max"] <- 0
    cvg_constraints(settings)$stop_codon["max"] <- 0
    # create various length primers
    constraints(settings)$primer_length <- c("min" = 18, "max" = 25) # was 18 to 28 for first 30 bases ... but reduced to preserve memory for larger region
    # create several optimal sets allowing for different penalty deviations:
    cur.results.loc <- file.path("results", paste0(sample.name, "_exploration"), analysis.type)
    dir.create(cur.results.loc)
    sample.name <- paste0(sample.name, "_", analysis.type)
    primer.data <- design_primers_global_single(template.df, settings, mode.directionality = mode.directionality,
                    sample.name = sample.name, target.temps = target.temps, init.algo = init.algo, 
                    max.degen = max.degen, conservation = conservation, cur.results.loc = cur.results.loc, prefilter = prefilter)
    template.data <- replicate(length(primer.data), template.df, simplify = FALSE)
    save(primer.data, file = file.path(cur.results.loc, paste0("primer_data_", sample.name, ".Rdata")))
    source("src/analyses/primer_comparison.R") # load comparison plot function
    create_comparison_plots(primer.data, template.data, settings, 
                            cur.results.loc, sample = "IGH", plot.details = FALSE)
    # select only 100% coverage sets
    cvg <- sapply(seq_along(primer.data), function(x) get_cvg_ratio(primer.data[[x]], template.data[[x]]))
    idx <- which(cvg == 1) # cvg == 1
    if (length(idx) == 0) {
        warning("No primer set with 100% coverage!")
    }
    create_comparison_plots(primer.data[idx], template.data[idx], settings, 
                            cur.results.loc, sample = "IGH_sel", plot.details = FALSE)

}

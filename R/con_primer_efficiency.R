##################
# Primer efficiency
###################
#' Primer Efficiency.
#'
#' Computes the efficiency of primer binding events for Taq polymerase.
#'
#' This function uses DECIPHER's \link[DECIPHER]{CalculateEfficiencyPCR}.
#'
#' @param fw.primers Primer sequence strings.
#' @param fw.start Binding position (start).
#' @param fw.end Binding position (end).
#' @param covered List of covered template indices per primer.
#' @param taqEfficiency Whether the efficiency shall be computed
#' using a mismatch-model developed for Taq polymerases. The default setting
#' is \code{TRUE}. Set \code{taqEfficiency} to \code{FALSE} if you are using
#' another polymerase than Taq.
#' @param annealing.temp Annealing temperature for which to evaluate efficiency.
#' @param primer_conc Primer concentration.
#' @param sodium.eq.concentration The sodium-equivalent concentration of ions.
#' @param mode.directionality Primer directionality.
#' @param seqs Template sequence strings.
#' @param eff.only Compute only 3' terminal efficiencies (no thermodynamic model).
#' @return The efficiencies of primer binding events.
#' @keywords internal
#' @references Wright, Erik S., et al. 
#' "Exploiting extension bias in polymerase chain reaction to improve 
#' primer specificity in ensembles of nearly identical DNA templates." 
#' Environmental microbiology 16.5 (2014): 1354-1365.
# note: since we compute 'efficiency' for individual primers only, this means that the model for primer-primer formation does not do anything, except computing something relating to self dimers. this is ok, since openPrimeR offers its own methods for cross dimerization.
compute.efficiency <- function(fw.primers, fw.start, fw.end, covered, 
                               taqEfficiency, annealing.temp, primer_conc, 
                               sodium.eq.concentration, 
                               mode.directionality = c("fw", "rev"), seqs, 
                               eff.only = FALSE) {
    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    ori.primers <- fw.primers
    # select only viable sequences for efficiency considerations
    sel.idx <- which(fw.primers != "")
    fw.primers <- my.disambiguate(DNAStringSet(fw.primers[sel.idx]))
    fw.effs <- vector("list", length(fw.primers))
    # parallelize along primers:
    i <- NULL
    #for (i in seq_along(fw.primers)) { # for debug ..
    fw.effs <- foreach (i = seq_along(fw.primers), .combine = c) %dopar% {
        #print(i)
        Ta <- annealing.temp[i]
        starts <- fw.start[[sel.idx[i]]]
        ends <- fw.end[[sel.idx[i]]]
        covered.seq.idx <- covered[[sel.idx[i]]]
        covered.seqs <- seqs[covered.seq.idx]
        cur.data <- data.frame(Primers = rep(NA, length(covered.seqs)), Templates = rep(NA, length(covered.seqs)))
        for (j in seq_along(starts)) {
            # identify template regions that are covered
            targets <- Biostrings::extractAt(covered.seqs[[j]], at = IRanges(start = starts[j], 
                end = ends[j]))  # fw ambig target sequences
            combis <- expand.grid(as.character(unlist(targets)), as.character(fw.primers[[i]]), stringsAsFactors = FALSE)  # possible binding modes for ambig sequences
            # select only the best (i.e. the actual) binding mode -> the one with smallest distance
            d <- stringdist::stringdist(combis[,1], combis[,2], nthread = 1)
            conformation.idx <- which.min(d)
            cur.data$Primers[j] <- combis[conformation.idx, 2]
            cur.data$Templates[j] <- combis[conformation.idx, 1]
        }
        # select only unique combinations of primers and templates
        dup.data <- unlist(lapply(seq_len(nrow(cur.data)), function(x) paste(unlist(cur.data[x,]), collapse = "")))
        # create an index to compute efficiency only for unique combinations of primers&templates
        idx <- which(!duplicated(dup.data))
        m <- match(dup.data, dup.data[idx]) # idx of corresponding representatitives
        if (nrow(cur.data) == 0) { # no templates covered
            cur.eff <- 0 # primer has 0 efficiency
        } else if (mode.directionality == "fw") {
            # reverse complement the binding sequence
            if (eff.only) {
                cur.eff <- .Call("terminalMismatch", cur.data$Primers[idx], rev.comp.sequence(cur.data$Templates[idx]), maxDistance = 1, 
                    maxGaps = 0, processors = 1, # set processors to 1 due to our external foreach loop!
                    PACKAGE = "openPrimeR")
            } else {
                cur.eff <- DECIPHER::CalculateEfficiencyPCR(DNAStringSet(cur.data$Primers[idx]), 
                        Biostrings::reverseComplement(
                            DNAStringSet(cur.data$Templates[idx])
                        ), Ta, primer_conc, 
                        sodium.eq.concentration, 
                        taqEfficiency = taqEfficiency, processors = 1)
            }
        } else if (mode.directionality == "rev") {
            if (eff.only) {
                cur.eff <- .Call("terminalMismatch", cur.data$Primers[idx], cur.data$Templates[idx], maxDistance = 1, 
                    maxGaps = 0, processors = 1,
                    PACKAGE = "openPrimeR")
            } else {
                # no reverse complement is necessary
                cur.eff <- DECIPHER::CalculateEfficiencyPCR(DNAStringSet(cur.data$Primers[idx]), 
                            DNAStringSet(cur.data$Templates[idx]), Ta, 
                            primer_conc, sodium.eq.concentration, 
                            taqEfficiency = taqEfficiency, processors = 1)
           }
        }
        eff <- rep(NA, nrow(cur.data))
        eff <- cur.eff[m]
        list(eff)
    }
    out.effs <- vector("list", length(ori.primers))
    out.effs[sel.idx] <- fw.effs
    return(out.effs)
}
#' Primer Efficiency.
#'
#' Computes the efficiency of primer binding events for Taq polymerase.
#'
#' This function uses DECIPHER's \code{\link[DECIPHER]{CalculateEfficiencyPCR}}.
#'
#' @param primer.df Primer data frame.
#' @param template.df Template data frame.
#' @param annealing.temp Annealing temperature for which to evaluate efficiency.
#' @param taqEfficiency Whether the efficiency shall be computed
#' using a mismatch-model developed for Taq polymerases. The default setting
#' is \code{TRUE}. Set \code{taqEfficiency} to \code{FALSE} if you are using
#' another polymerase than Taq.
#' @param primer_conc Primer concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param eff.only Compute only 3' terminal efficiencies (no thermodynamic model).
#' @param mode Compute efficiencies for on-target coverage events (\code{on_target})
#' or off-target coverage events (\code{off_target}).
#' @return A list with the efficiency of every primer binding event.
#' @keywords internal
#' @examples
#' data(Ippolito)
#' p <- PCR(settings)
#' # Requires OligoArrayAux software:
#' \dontrun{
#' eff.df <- compute.primer.efficiencies(primer.df, template.df, 55, 
#'              p$primer_concentration, p$Na_concentration,
#'              p$Mg_concentration, p$K_concentration, p$Tris_concentration)
#' } 
compute.primer.efficiencies <- function(primer.df, template.df, annealing.temp, 
    taqEfficiency, primer_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, 
    eff.only = FALSE, mode = c("on_target", "off_target")) {

    if (!"primer_coverage" %in% colnames(primer.df)) {
        stop("Efficiency computations require primer coverage.")
    }
    if (!check.tool.function()["OligoArrayAux"]) {
        # cannot test without OligoArrayAux
        stop("Cannot compute efficiency: OligoArrayAux not available.")
    }
    if (length(which(!is.na(annealing.temp))) != nrow(primer.df)) {
        stop("Annealing temperatures not provided for all primers.")
    }
    #message("Computing efficiency @ annealing temp: ", annealing.temp)
    #time <- Sys.time()
    if (mode == "on_target") {
        fw.start <- lapply(primer.df$Binding_Position_Start_fw, function(x) as.numeric(unlist(strsplit(x, 
            split = ","))))
        fw.end <- lapply(primer.df$Binding_Position_End_fw, function(x) as.numeric(unlist(strsplit(x, 
            split = ","))))
        rev.start <- lapply(primer.df$Binding_Position_Start_rev, function(x) as.numeric(unlist(strsplit(x, 
            split = ","))))
        rev.end <- lapply(primer.df$Binding_Position_End_rev, function(x) as.numeric(unlist(strsplit(x, 
            split = ","))))
        covered <- lapply(primer.df$Covered_Seqs, function(x) match(as.numeric(unlist(strsplit(x, 
            split = ","))), template.df$Identifier))
    } else {
        fw.start <- lapply(primer.df$Off_Binding_Position_Start_fw, function(x) as.numeric(unlist(strsplit(x, 
            split = ","))))
        fw.end <- lapply(primer.df$Off_Binding_Position_End_fw, function(x) as.numeric(unlist(strsplit(x, 
            split = ","))))
        rev.start <- lapply(primer.df$Off_Binding_Position_Start_rev, function(x) as.numeric(unlist(strsplit(x, 
            split = ","))))
        rev.end <- lapply(primer.df$Off_Binding_Position_End_rev, function(x) as.numeric(unlist(strsplit(x, 
            split = ","))))
        covered <- lapply(primer.df$Off_Covered_Seqs, function(x) match(as.numeric(unlist(strsplit(x, 
            split = ","))), template.df$Identifier))
    }
    seqs <- my.disambiguate(DNAStringSet(template.df$Sequence))
    sodium.eq.concentration <- compute.sodium.equivalent.conc(na_salt_conc, 
                               mg_salt_conc, k_salt_conc, tris_salt_conc)
    fw.effs <- compute.efficiency(primer.df$Forward, fw.start, fw.end, covered, 
                                  taqEfficiency, annealing.temp, primer_conc, 
                                  sodium.eq.concentration,"fw", seqs, eff.only = eff.only)
    rev.effs <- compute.efficiency(primer.df$Reverse, rev.start, rev.end, 
                                   covered, taqEfficiency, annealing.temp, primer_conc, 
                                   sodium.eq.concentration, "rev", seqs, eff.only = eff.only)
    #message("Runtime: ", Sys.time() - time) 
    # combine efficiciencies from fw and rev primers if necessary
    fw.idx <- which(primer.df$Forward != "")
    rev.idx <- which(primer.df$Reverse != "")
    effs <- lapply(seq_len(nrow(primer.df)), function(j) {
        if (j %in% fw.idx && j %in% rev.idx) {
            # compute the geometric mean of forward and reverse efficiencies
            sqrt(fw.effs[[j]] * rev.effs[[j]])
        } else if (j %in% fw.idx) {
            fw.effs[[j]]
        } else if (j %in% rev.idx) {
            rev.effs[[j]]
        }
    })
    return(effs)
}

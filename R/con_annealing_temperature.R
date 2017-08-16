######
# Constraint: annealing temperature
########
#' Annealing temperature.
#'
#' Identifies the optimal annealing temperature of a set of primers.
#' If primers cover template sequences, the annealing temperature
#' is computed using Rychlik's formula.
#' Otherwise, the annealing temperature is determined using the
#' rule of thumb based on the melting temperatures of the primers.
#'
#' @param primer.df Primer data frame.
#' @param mode.directionality Primer directionality.
#' @param template.df Template data frame
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param primer_conc Primer concentration.
#' @return The optimal annealing temperature.
#' @keywords internal
compute_annealing_temp <- function(primer.df, mode.directionality, 
                      template.df, na_salt_conc, mg_salt_conc, k_salt_conc, 
                      tris_salt_conc, primer_conc) {
    Ta <- NULL
    #if ("primer_coverage" %in% colnames(primer.df)) {
        ## "covered Ta" only makes sense for non-mismatch binding events ...? quite low values ...
        #Ta <- compute.covered.Ta(primer.df, mode.directionality, 
                          #template.df, na_salt_conc, mg_salt_conc, k_salt_conc, 
                          #tris_salt_conc, primer_conc)
        ## take the smallest annealing temperature
        #if (length(Ta) != 0) {
            #Ta <- min(Ta, na.rm = TRUE)
        #}
    #}
    #if (length(Ta) == 0) { # no primer cvg
        # -> use the rule of thumb based on Tm
    #} 
    compute.idx <- seq_len(nrow(primer.df))
    if ("melting_temp" %in% colnames(primer.df)) {
        compute.idx <- which(is.na(primer.df$melting_temp))
    }
    if (length(compute.idx) != 0) {
        # compute melting temps
         melting.temps <- compute.melting.temps(primer.df[compute.idx,], primer_conc, 
            na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, 
            mode.directionality)
        p.df <- update.constraint.values(primer.df[compute.idx, ], melting.temps)
        if ("melting_temp" %in% colnames(primer.df)) {
            # overwrite existing MELT entries with new rows
            primer.df[compute.idx, ] <- p.df
        } else {
            # add the columns to the data frame
            primer.df <- cbind(primer.df, p.df)
        }
    } 
    Ta <- annealing.temp.rule.of.thumb(primer.df$melting_temp)
    return(Ta)
}

#' Rule of thumb for annealing temperature
#'
#' Computes the annealing temperature using a rule of thumb
#' 
#' @param melting.temp Melting temperatures of primers
#' 
#' @return The annealing temperature corresponding to the input melting temperature.
#' @keywords internal
annealing.temp.rule.of.thumb <- function(melting.temp) {
    if (length(melting.temp) == 0 || all(is.na(melting.temp ))) {
        return(NA)
    }
    Ta <- melting.temp - 5
    return(Ta)
}

#' Annealing temperature
#'
#' Computes the annealing temperature using all binding events.
#'
#' @param primer.df Primer data frame.
#' @param mode.directionality Primer directionality.
#' @param template.df Template data frame
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param primer_conc Primer concentration.
#'
#' @return The recommended annealing temperature.
#' @keywords internal
compute.covered.Ta <- function(primer.df, mode.directionality = c("fw", "rev", "both"), template.df, na_salt_conc, 
                               mg_salt_conc, k_salt_conc, tris_salt_conc, primer_conc) {
    if (length(mode.directionality) == 0) {
        stop("'mode.directionality' not supplied.")
    }
    warning("covered.Ta is deprecated: Ta is too low!")
    mode.directionality <- match.arg(mode.directionality)
    if (!"Covered_Seqs" %in% colnames(primer.df) || length(primer.df) == 0) {
        stop("Cannot compute annealing temperature: primer coverage is not available yet!")
    }
    seq.idx <- covered.seqs.to.idx(primer.df$Covered_Seqs, template.df)
    Ta <- na.omit(unlist(lapply(1:nrow(primer.df), function(x) {
        compute.Ta(primer.df[x, ], template.df[seq.idx[[x]], ], 
                   mode.directionality, na_salt_conc, mg_salt_conc, 
                   k_salt_conc, tris_salt_conc, primer_conc)
        })))
    return(Ta)
}
#' Annealing temperature
#'
#' Computes the annealing temperature using all binding events.
#'
#' @param primer.df Primer data frame.
#' @param template.df Template data frame
#' @param mode.directionality Primer directionality.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param primer_conc Primer concentration.
#'
#' @return All annealing temperatures for given binding events.
#' @references 
#' Rychlik, W. J. S. W., W. J. Spencer, and R. E. Rhoads. "Optimization of the annealing temperature for DNA amplification in vitro." Nucleic acids research 18.21 (1990): 6409-6412.
#' @keywords internal
compute.Ta <- function(primer.df, template.df, mode.directionality = c("fw", "rev", "both"), 
                       na_salt_conc, mg_salt_conc, 
                       k_salt_conc, tris_salt_conc, primer_conc) {
    if (length(mode.directionality) == 0) {
        stop("'mode.directionality' not supplied.")
    }
    mode.directionality <- match.arg(mode.directionality)
    if (length(template.df) == 0 || nrow(template.df) == 0) {
        return(NA)
    }
    if (!"melting_temp" %in% colnames(primer.df)) {
        Tm.p <- compute.melting.temps(primer.df, primer_conc, na_salt_conc, mg_salt_conc, 
                k_salt_conc, tris_salt_conc, mode.directionality)$melting_temp  # primer Tm
    } else {
        Tm.p <- primer.df$melting_temp
    }
    if (!"melting_temp" %in% colnames(template.df)) {
        Tm.s <- compute.Tm.baldino(template.df$Sequence, na_salt_conc, mg_salt_conc, k_salt_conc, 
            tris_salt_conc, primer_conc)  # template Tm
    } else {
        Tm.s <- template.df$melting_temp
    }
    result <- 0.3 * Tm.p + 0.7 * Tm.s - 14.9
    na.idx <- which(is.nan(result))
    if (length(na.idx) != 0) {
        result[na.idx] <- NA
    }
    return(result)
}


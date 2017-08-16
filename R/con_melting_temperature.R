#######
# Melting temperature computations
######

#' Change in Free Energy.
#'
#' Computes the change in free energy.
#'
#' @param delta.H Change in enthalpy in cal/mol.
#' @param delta.S Change in entropy cal/mol*K.
#' @param temp Temperature in Celsius for which to compute free energy change.
#'
#' @return The change in free energy in kcal/mol.
#' @keywords internal
get.delta.G <- function(delta.H, delta.S, temp = 37) {
    K <- convert.temperature(temp, "K")
    G <- (delta.H - (K * delta.S))/1000
    return(G)  # in kcal/mol
}
#' Conversion from J to cal
#'
#' Converts the input from Joule to calories.
#'
#' @param val.J Numeric Joule value.
#'
#' @return The value correspdoning to \code{val.J} in calories.
#' @keywords internal
joule.to.cal <- function(val.J) {
    return(val.J/4.184)
}
#' Conversion between Celsius and Kelvin
#'
#' Converts the input from Kelvin to Celsius or from Celsius to Kelvin.
#'
#' If \code{temp.scale} is 'K', \code{T_m} is transformed from Celsius
#' to Kelvin. If \code{temp.scale} is 'C', \code{T_m} is transformed from
#' Kelvin to Celsius. The default is to transform from Celsius to Kelvin.

#' @param temp The input temperature.
#' @param temp.scale The desired unit of the output temperature.
#'
#' @return Transforms the input temperature to the specified \code{temp.scale}.
#' @keywords internal
convert.temperature <- function(temp, temp.scale = c("K", "C")) {
    # scale: kelvin -> convert to kelvin scale: celsius -> convert to celsius
    if (length(temp.scale) == 0) {
        stop("Please supply the 'temp.scale' parameter.")
    }
    temp.scale <- match.arg(temp.scale)
    if (temp.scale == "K") {
        return(temp + 273.15)
    } else if (temp.scale == "C") {
        return(temp - 273.15)
    } 
}

#' Baldino Formula
#'
#' Computes the melting temperature using the formulation by Baldino.
#'
#' @param sequences Input sequence strings.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param primer_conc Primer concentration.
#'
#' @return The melting temperature for the input sequences.
#' @references 
#' Rychlik, W. J. S. W., W. J. Spencer, and R. E. Rhoads. "Optimization of the annealing temperature for DNA amplification in vitro." Nucleic acids research 18.21 (1990): 6409-6412.
#' @keywords internal
compute.Tm.baldino <- function(sequences, na_salt_conc, mg_salt_conc, k_salt_conc, 
    tris_salt_conc, primer_conc) {
    l <- nchar(sequences)  # assume longest length of the PCR product
    gc.ratio <- compute.gc.ratio(sequences)
    sodium.eq.concentration <- compute.sodium.equivalent.conc(na_salt_conc, mg_salt_conc, 
                                            k_salt_conc, tris_salt_conc)  
    T_mC <- (0.41 * (gc.ratio * 100)) + 
            (16.6 * log10(sodium.eq.concentration*1000)) - (675/l)
    return(T_mC)
}
#' Non-Thermodynamic Computation of Melting Temperatures.
#'
#' Computes the melting temperature of primers from an empiric formula.
#'
#' @param primer.df A \code{Primers} object.
#' @return A data frame with melting temperature information for the primers.
#' @keywords internal
compute.empiric.melting.temp <- function(primer.df) {
    fw.seqs <- convert.from.iupac(primer.df$Forward)
    rev.seqs <- convert.from.iupac(primer.df$Reverse)
    Tm <- function(x) {
        if (x[[1]] == "") {
            return(NA)
        }
        x <- strsplit(x, split = "")
        A <- sapply(x, function(y) length(which(y == "a")))
        C <- sapply(x, function(y) length(which(y == "c")))
        G <- sapply(x, function(y) length(which(y == "g")))
        T <- sapply(x, function(y) length(which(y == "t")))
        # return the smallest reuired melting temperature
        Tm <- 64.9 + 41 * (G + C - 16.4) / (A + T + G + C) 
        return(min(Tm))
    }
    Tm.fw <- sapply(fw.seqs, Tm)
    Tm.rev <- sapply(rev.seqs, Tm)
    Tm.K.fw <- convert.temperature(Tm.fw, "K")
    Tm.K.rev <- convert.temperature(Tm.rev, "K")
    # return overall minimal Tm of both directions:
    melting.temp <- sapply(seq_len(nrow(primer.df)), function(x) min(Tm.fw[x], Tm.rev[x], na.rm = TRUE))
    Tm.data <- data.frame("melting_temp" = melting.temp, "Tm_C_fw" = Tm.fw, "Tm_C_rev" = Tm.rev, "Tm_K_fw" = Tm.K.fw, "Tm_K_rev" = Tm.K.rev)
    return(Tm.data)
}

#' Computation of Thermodynamic Melting Temperatures.
#'
#' Use nearest-neighbor thermodynamic computations to find melting temperatures.
#'
#' @param primer.df Primer data frame.
#' @param primer_conc Primer concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param mode.directionality Direction of primers
#' @return Data frame with melting temperature info for the input primers.
#' @keywords internal
compute.melting.temps.thermo <- function(primer.df, primer_conc, na_salt_conc, mg_salt_conc, 
    k_salt_conc, tris_salt_conc, mode.directionality = c("fw", "rev", "both")) {

    used.Tm <- rep(NA, nrow(primer.df))  # Tm to be used in optimization algorithm
    Tm.fw <- call.melt(primer.df$Forward, complement.sequence(primer.df$Forward), 
        primer_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc)[, 3:6]
    Tm.rev <- call.melt(primer.df$Reverse, complement.sequence(primer.df$Reverse), 
        primer_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc)[, 3:6]
    cnames <- NULL
    if (mode.directionality == "fw" || mode.directionality == "both" && length(Tm.fw) != 0) {
        cnames <- colnames(Tm.fw)
    } else {
        cnames <- colnames(Tm.rev)
    }
    # set empty data frame to NA
    if (length(Tm.fw) != 0) {
        colnames(Tm.fw) <- paste(cnames, "_fw", sep = "")
    } else if (length(Tm.fw) == 0 && length(Tm.rev) != 0) {
        Tm.fw <- Tm.rev
        Tm.fw[TRUE] <- NA
        colnames(Tm.fw) <- paste(cnames, "_fw", sep = "")
    }
    # set empty data frame to NA
    if (length(Tm.rev) != 0) {
        colnames(Tm.rev) <- paste(cnames, "_rev", sep = "")
    } else if (length(Tm.rev) == 0 && length(Tm.fw) != 0) {
        Tm.rev <- Tm.fw
        Tm.rev[TRUE] <- NA
        colnames(Tm.rev) <- paste(cnames, "_rev", sep = "")
    }
    if (length(Tm.fw) != 0 && length(Tm.rev) != 0) {
        used.Tm <- unlist(lapply(1:nrow(Tm.fw), function(x) min(c(Tm.fw$Tm_C_fw[x], Tm.rev$Tm_C_rev[x]), 
            na.rm = TRUE)))  # use the min Tm 
    } else if (length(Tm.fw) != 0 && length(Tm.rev) == 0) {
        used.Tm <- Tm.fw$Tm_C_fw
    } else if (length(Tm.rev) != 0 && length(Tm.fw) == 0) {
        used.Tm <- Tm.rev$Tm_C_rev
    } else {
        used.Tm <- NULL
    }
    Tm.data <- cbind(Tm.fw, Tm.rev, melting_temp = used.Tm)
    return(Tm.data)
}
#' Computation of Melting Temperatures.
#'
#' Use nearest-neighbor thermodynamic computations to find melting temperatures.
#'
#' @param primer.df Primer data frame.
#' @param primer_conc Primer concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param mode.directionality Direction of primers
#' @return Data frame with melting temperature info for the input primers.
#' @keywords internal
compute.melting.temps <- function(primer.df, primer_conc, na_salt_conc, mg_salt_conc, 
    k_salt_conc, tris_salt_conc, mode.directionality = c("fw", "rev", "both")) {

    if (length(mode.directionality) == 0) {
        stop("Please supply the 'mode.directionality' argument.")
    }
    mode.directionality <- match.arg(mode.directionality)
    if (check.tool.function()["MELTING"]) {
        # melting is available -> use thermodynamic model
        Tm.data <- compute.melting.temps.thermo(primer.df, primer_conc, na_salt_conc, mg_salt_conc, 
                            k_salt_conc, tris_salt_conc, mode.directionality)
    } else {
        # melting is not available -> use empiric formula without salt correction
        Tm.data <- compute.empiric.melting.temp(primer.df)
    }
    # determine max temperature diff for every primer
    temp.diff <- get.melting.temp.diff(Tm.data$Tm_C_fw, Tm.data$Tm_C_rev)
    result <- cbind(Tm.data,melting_temp_diff = temp.diff)
    return(result)
}
#' Computation of Maximal Melting Temperature Differences.
#'
#' @param Tm.fw The melting temperatures of forward primers.
#' @param Tmr.rev The melting temperatures of reverse primers.
#' @return The worst-case melting temperature difference, for every primer.
#' @keywords internal
get.melting.temp.diff <- function(Tm.fw, Tm.rev) {
    if (length(Tm.fw) != 0) {
        fw.idx <- which(!is.na(Tm.fw))
    } else {
        fw.idx <- NULL
    }
    if (length(Tm.rev) != 0) {
        rev.idx <- which(!is.na(Tm.rev))
    } else {
        rev.idx <- NULL
    }
    max.len <- max(length(Tm.fw), length(Tm.rev))
    temp.diff.fw <- rep(0, max.len)
    temp.diff.rev <- rep(0, max.len)
    if (length(fw.idx) != 0) {
        temp.diff.fw[fw.idx] <- sapply(fw.idx, function(x) max(c(abs(Tm.fw[x] - Tm.fw), abs(Tm.fw[x] - Tm.rev)), na.rm = TRUE))
    }
    if (length(rev.idx) != 0) {
        temp.diff.rev[rev.idx] <- sapply(rev.idx, function(x) max(c(abs(Tm.rev[x] - Tm.rev), abs(Tm.rev[x] - Tm.fw)), na.rm = TRUE))
    }
    temp.diff <- sapply(seq_along(temp.diff.fw), function(x) max(c(temp.diff.fw[x], temp.diff.rev[x])))
    return(temp.diff)
}
#' Thermodynamic melting temperature computations.
#'
#' Computes the melting temperature for the input primers.
#'
#' @param primers List of primer strings.
#' @param complements List with corresponding complements.
#' @param out.file Path to the file where MELTING will write the results.
#' @param primer_conc Primer concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @param ID identifiers of the input primers
#' @return Melting temperature data frame. 
#' @keywords internal
call.melt.single <- function(primers, complements, out.file, primer_conc, na_salt_conc, 
    mg_salt_conc, k_salt_conc, tris_salt_conc, ID) {
    # call melt command ID: used to keep track of ambiguous seqs that were generated
    if (length(primers) == 0 || length(complements) == 0) {
        return(NULL)
    }
    melt <- Sys.which("melting-batch")
    if (melt == "") {
        stop("MELTING not available.")
    } else {
		if (grepl("windows", .Platform$OS.type)) {
			# for windows, use the bat file
			jar.file <- Sys.which("melting5.jar")
			melt <- paste("java -cp", jar.file, "melting.BatchMain")
		}
		#message(paste("Using the following melting location: ", melt))
	}
    # melt results from b atch are in joul -> convert to cal use Allawi, SantaLucia
    # 1997 publication nearest neighbor data
    nn <- "all97"  # nearest neighbor model
    am <- "wetdna91"  # from Wetmur1991, used for more than 60 bases primers
    # standard parameters: -P: nucleic acid strand in excess (template or primers) -E
    # <compound>: concentration of compound in molar if [Mg] is 0, don't add the
    # parameter
    salt.cmd <- paste(" -E Na=", format(na_salt_conc, scientific = FALSE), sep = "")
    if (mg_salt_conc != 0) {
        salt.cmd <- paste(salt.cmd, " -E Mg=", format(mg_salt_conc, scientific = FALSE), 
            sep = "")
    }
    if (k_salt_conc != 0) {
        salt.cmd <- paste(salt.cmd, " -E K=", format(k_salt_conc, scientific = FALSE), 
            sep = "")
    }
    if (tris_salt_conc != 0) {
        salt.cmd <- paste(salt.cmd, " -E Tris=", format(tris_salt_conc, scientific = FALSE), 
            sep = "")
    }
    # don't use custom tandem file anymore, not necessary if we don't do mismatch Tm computation.
    #tandem.mm.file <- "AllawiSantaluciaPeyret1997_1998_1999tanmm_mod.xml"
    melt.cmd <- paste(melt, " -am ", am, " -nn ", nn, " -H dnadna -P ", format(primer_conc, 
        scientific = FALSE), salt.cmd, #" -tan :", tandem.mm.file, 
        " ", out.file, sep = "")
    #message(melt.cmd)
    system(melt.cmd, ignore.stdout = TRUE)  # don't igore stdout for debugging TODO, set to TRUE for debugging
    # read results
    result.file <- paste(out.file, ".results.csv", sep = "")
    raw.data <- readChar(result.file, file.info(result.file)$size)
    # modify the header change a separator from \r to \t
    raw.data <- gsub("\r", "\t", raw.data)
	# remove terminal tabulators:
	raw.data <- gsub("\t\n", "\n", raw.data)
    results <- try(read.delim(text = raw.data, header = TRUE, 
                stringsAsFactors = FALSE), silent = FALSE)
    # raw data is empty sometime. why if there's input -> melting fails i guess?
    if (class(results) == "try-error") {
        warning("MELTING could not compute the Tm:\n", 
            "File: ", result.file, "\n", 
            attr(results, "condition"))
        e <- rep(NA, length(primers))
        return(data.frame(Sequence = e, Complementary = e, DeltaH = e, DeltaS = e, 
            Tm_C = e, Tm_K = e, ID = ID))
    }
    # determine system locale settings:
    decimal.sep <- Sys.localeconv()["mon_decimal_point"]
    if (decimal.sep == ".") {
        decimal.sep <- "\\."
        mill.sep <- ","
    } else {
        mill.sep <- "\\."
    }
    # replace 1000pos commas and convert from J to calories
    results[, 3] <- joule.to.cal(as.numeric(gsub(decimal.sep, ".", gsub(mill.sep, 
        "", results[, 3]))))
    results[, 4] <- joule.to.cal(as.numeric(gsub(decimal.sep, ".", gsub(mill.sep, 
        "", results[, 4]))))
    results[, 5] <- as.numeric(gsub(decimal.sep, ".", gsub(mill.sep, "", results[, 
        5])))
    colnames(results)[5] <- "Tm_C"
    results <- cbind(results, Tm_K = convert.temperature(results$Tm_C, "K"))
    ######### insert missing sequence results, e.g. for '' seqs
    sel <- rep(NA, length(primers))
    for (i in seq_along(primers)) {
        primer <- primers[i]
        comp <- complements[i]
        idx.s <- which(sapply(results$Sequence, function(x) x == primer))
        idx.c <- which(sapply(results$Complementary, function(x) x == comp))
        idx <- intersect(idx.s, idx.c)
        if (length(idx) != 0) {
            sel[i] <- idx[1]
        }
    }
    results <- results[sel, ]
    # insert sequence and complementary for missing entries
    missing.idx <- which(is.na(results$Sequence))
    if (length(missing.idx) != 0) {
        results[missing.idx, "Sequence"] <- primers[missing.idx]
        results[missing.idx, "Complementary"] <- complements[missing.idx]
    }
	# prepare output:
    results <- cbind(results, ID = ID)
    results$Sequence <- tolower(results$Sequence)
    results$Complementary <- tolower(results$Complementary)
    return(results)
}

#' Thermodynamic melting temperature computations.
#'
#' Computes the melting temperature for the input primers.
#'
#' @param primers Character vector of primer strings.
#' @param complements Character vector with complement sequences
#' corresponding to \code{primers}.
#' @param primer_conc Primer concentration.
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @return Melting temperature data frame. 
#' @keywords internal
#' @references Le Novère N. (2001). MELTING, computing the 
#' melting temperature of nucleic acid duplex. Bioinformatics, 17: 1226-1227. 
#' Dumousseau M., Rodriguez N., Juty N., Le Novère N. (2012) MELTING, 
#' a flexible platform to predict the melting temperatures of nucleic acids. 
#' BMC Bioinformatics, 13: 101.
call.melt <- function(primers, complements, primer_conc, na_salt_conc, mg_salt_conc, 
    k_salt_conc, tris_salt_conc) {
    # complement needs to be specified for mismatch consideration 
    # ./melting -S TGAGGTGCAGCTGGTGGAGTC -H dnadna -P 0.00000005 -E Na=0.05 -v -C ACTCCACGCCGACCACCTCAA 
    # ./melting-batch -H dnadna -P 0.00000005 -E Na=0.05
    # primer_rev_comp.fasta complement from 3' to 5' 1. create file for melt
    if (all(primers == "") || length(primers) == 0 || length(complements) == 0) {
        return(NULL)
    }
    #print("call.melt: primer input")
    #print(primers)
    runID <- digest::digest(c(primers, complements), "md5")
    if (length(primers) != length(complements)) {
        stop("The number of primers did not correspond to the number of complements.")
    }
    if (Sys.which("melting-batch") == "") {
        stop("MELTING executable not available. Please install MELTING to compute melting temperatures.")
    }
    idx <- which(primers != "")
    out.len <- length(primers)
    primers <- primers[idx]
    complements <- complements[idx]
    # disambiguate primers
    unambig.primers <- convert.from.iupac(primers)
    # determine idx of exploded primers
    unambig.primers.c <- convert.from.iupac(complements)  # unambiguated complements
    # need to get all combinations of unambig primers with unambig.primers.c
    ambig.combinations <- lapply(seq_along(unambig.primers), function(x) expand.grid(unambig.primers[[x]], 
        unambig.primers.c[[x]], stringsAsFactors = FALSE))
    ids <- unlist(lapply(seq_along(ambig.combinations), function(x) rep(x, nrow(ambig.combinations[[x]]))))
    out.primers <- unlist(lapply(seq_along(ambig.combinations), function(x) toupper(ambig.combinations[[x]][, 
        1])))
    out.complements <- unlist(lapply(seq_along(ambig.combinations), function(x) toupper(ambig.combinations[[x]][, 
        2])))
    # MELTING can't handle terminal mismatches, remove those
    check.comp <- toupper(complement.sequence(out.complements))
    rm.idx <- which(substring(check.comp, 1, 1) != substring(out.primers, 1, 1) |
        substring(check.comp, nchar(check.comp), nchar(check.comp)) !=
        substring(out.primers, nchar(out.primers), nchar(out.primers)))
    if (length(rm.idx) != 0) {
        out.primers <- out.primers[-rm.idx]
        ids <- ids[-rm.idx]
        out.complements <- out.complements[-rm.idx]
    }
    melt.tab <- data.frame(Identifier = ids, Primer = out.primers, Complement = out.complements)
    # write one file for each parallel call
    seqs.per.file <- ceiling(length(out.primers)/foreach::getDoParWorkers())
    no.files <- ceiling(length(out.primers)/seqs.per.file)
    out.files <- rep(NA, no.files)
    total.seqs.covered <- 0
    seq.idx <- vector("list", no.files)  # idx of the primers in out.primers per file
    # write out.files
    for (i in seq_len(no.files)) {
        seqs.added <- min(seqs.per.file, length(out.primers) - total.seqs.covered)
        s <- NULL
        if (i == 1) {
            s <- 1
        } else {
            s <- max(seq.idx[[i - 1]]) + 1
        }
        e <- s + (seqs.added - 1)
        seq.idx[[i]] <- s:e
        total.seqs.covered <- total.seqs.covered + seqs.added
        out.file <- tempfile(paste("melt_seqs_", i, runID, sep = ""), fileext = ".txt")
        out.files[i] <- out.file
        #print("melt Seqs:")
        #print(melt.tab[seq.idx[[i]],])
        write.table(melt.tab[seq.idx[[i]], c("Primer", "Complement")], file = out.file, 
					sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
    }
    # parallelized call:
    i <- NULL
    results <- foreach(i = seq_along(out.files), .combine = "rbind") %dopar% {
        result <- call.melt.single(out.primers[seq.idx[[i]]], out.complements[seq.idx[[i]]], 
            out.files[i], primer_conc, na_salt_conc, mg_salt_conc, k_salt_conc, tris_salt_conc, 
            ids[seq.idx[[i]]])
    }
    # break results down to individual sequences using the ID column to integrate
    # ambig seq results
    u.IDs <- unique(results$ID)
    results$deltaG <- get.delta.G(results$DeltaH, results$DeltaS)
    # integrate by selecting worst-case deltaG per ambig seq combination
    results <- ddply(results, c("ID"), function(x) arrange(x, substitute(deltaG))[1, ])  # select unique pairs
    # replace Sequence/Complementary stretch with original input (ambiguous)
    results[, "Sequence"] <- primers[results$ID]
    results[, "Complementary"] <- complements[results$ID]
    results <- results[, colnames(results) != "ID"]  # remove ID column
    # create entries for primers without this direction ("")
    out <- results[1,]
    out[,] <- NA
    out <- do.call(rbind, replicate(out.len, out, simplify = FALSE))
    out[idx,] <- results
	# clean up temporary files
	on.exit(file.remove(out.files))
    return(out)
}

#' Sodium-equivalent Concentration
#'
#' Computes the sodium-equivalent concentration for the input ion concentrations.
#'
#' @param na_salt_conc Sodium ion concentration.
#' @param mg_salt_conc Magensium ion concentration.
#' @param k_salt_conc Potassium ion concentration.
#' @param tris_salt_conc Tris ion concentration.
#' @return The sodium-equivalent concentration of the input ion concentrations.
#' @references
#' Record, M. Thomas. "Effects of Na+ and Mg++ ions on the helix–coil transition of DNA." 
#' Biopolymers 14.10 (1975): 2137-2158.
#'
#' Owczarzy, Richard, et al. "Predicting stability of DNA duplexes in solutions containing magnesium and monovalent cations." 
#' Biochemistry 47.19 (2008): 5336-5353.
#'
#' Peyret, Nicolas. Prediction of nucleic acid hybridization: parameters and algorithms. 
#' Detroit: Wayne State University, 2000.
#' @keywords internal
compute.sodium.equivalent.conc <- function(na_salt_conc, mg_salt_conc, k_salt_conc, 
    tris_salt_conc) {
    # beta: correction factor from Peyret (2000)
    beta <- 3.3
    na_salt_eq <- beta * sqrt(mg_salt_conc) + na_salt_conc + k_salt_conc + tris_salt_conc/2  # sodium equivalent concentration
    return(na_salt_eq)
}


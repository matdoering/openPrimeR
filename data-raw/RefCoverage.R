##########
# Create the feature matrix for supervised learning using the experimentally determined coverage
#######

cellColor <- function(style) {
    fg  <- style$getFillForegroundXSSFColor()
    #fg <- style$getFillForegroundColor() # TODO wrong
    rgb <- tryCatch(fg$getRgb(), error = function(e) NULL)
    rgb <- paste(rgb, collapse = "")
    return(rgb)
}
get.ref.df <- function(color.df) {
    library(stringr)
    cnam <- colnames(color.df)
    counts <- str_count(cnam, "\\.")
    idx.1 <- which(counts == 1)
    idx.2 <- which(counts == 2)
    idx.3 <- which(counts == 3)
    cnam[idx.1] <- sub("\\.", "-", cnam[idx.1])
    cnam[idx.2] <- sub("\\.", "*", sub("\\.", "-", cnam[idx.2]))
    cnam[idx.3] <- sub("\\.", "*", sub("\\.", "-", sub("\\.", "-", cnam[idx.3])))
    m <- match(cnam, template.df$ID)
    # manual check:
    cbind(cnam, template.df$ID[m])
    color.df <- data.frame(color.df[,1], color.df[, !is.na(m)])
    colnames(color.df) <- c("Primer", cnam[!is.na(m)])
    return(color.df)
}
get.color.df <- function(xls.file, sheet.idx, nbr.rows) {
###############
# NEED TO RUN THIS ON WINDOWS!!!
###############
	library(xlsx)
	wb <- loadWorkbook(xls.file)
	sheet1 <- getSheets(wb)[[sheet.idx]] # first single primer evaluation, 18
	# get all rows
	rows  <- getRows(sheet1, seq_len(nbr.rows)) # if rows are not specified, too many are extracted
	cells <- getCells(rows)
	nbr.cols <- as.numeric(tail(strsplit(names(cells), split = "\\."), n =1)[[1]][2])
	styles <- vector("list", length(rows))
	for (i in seq_along(rows)) {
		cell.id <- paste0(i, ".", 1:nbr.cols)
		row.data <- cells[cell.id]
		cur.styles <- sapply(unlist(row.data), getCellStyle)
		styles[[i]] <- cur.styles
	}
	fill.colors <- lapply(styles, function(x) sapply(unlist(x), cellColor)) # if this is not the same length, add Identifier column leftmost
	fill.matrix <- do.call(rbind, fill.colors[2:length(fill.colors)]) # ignore header row
	col.df <- readColumns(sheet1, 1, nbr.cols, 1, nrow(fill.matrix) + 1) # doesn't read cell styles, nrow+1 to load header row
	cnames <- colnames(col.df)
	col.df[, 2:ncol(col.df)] <- fill.matrix[, 2:ncol(fill.matrix)] # first column has primer identifers -> don't overwrite!
	colnames(col.df) <- cnames
	return(col.df)
}
get_ref_data <- function(xls.file) {
    message("This function should be run on Windows for correct xls parsing ...")
    # turn color xls to labeled R data frame with colors indicating amplification (red: no, green: yes)
    tiller.ref.df <- get.color.df(xls.file, 1,13)
    tiller.ref.df <- get.ref.df(tiller.ref.df)
    open.ref.df <- get.color.df(xls.file, 2, 49)
    open.ref.df <- get.ref.df(open.ref.df)
    ref.data <- list(tiller.ref.df, open.ref.df)
    names(ref.data) <- c("Tiller2008", "openPrimeR")
    # compute tool cvg for tiller and openPrimeR
    primer.file.t <- system.file("extdata", "IMGT_data", "primers", 
                               "IGHV", "Tiller2008_1st.fasta",
                               package = "openPrimeR")
    primer.file.o <- system.file("extdata", "IMGT_data", "primers", 
                               "IGHV", "openPrimeR2017.fasta",
                               package = "openPrimeR")
    primer.files <- c(primer.file.t, primer.file.o)
    for (i in seq_along(primer.files)) {
        primer.df <- read_primers(primer.files[[i]])
        ref.data[[i]]$Primer <- unlist(lapply(primer.df$ID, function(x) rep(x, 3))) # annotate primers with right IDs (order MUST be the same as the provided 'primer.df' for correct cycling!) -> 3 because of triplicates
    }
    return(ref.data)
}
get_tool_primer_data <- function(primer.location, template.df, cur.settings) {
    # load primers
    primer.df <- openPrimeR::read_primers(primer.location, "_fw", "_rev")
    primer.df <- openPrimeR::check_constraints(primer.df, template.df, cur.settings)
    template.df <- openPrimeR::update_template_cvg(template.df, primer.df)
    out <- list("Primers" = primer.df, "Templates" = template.df)
    return(out)
}
get_learning_matrix <- function(ref.data) {
    data(Ippolito) # load templates
    # primer set to use for calibration:
    mode.directionality <- "fw"
    # create tool data for Tiller / openPrimeR:
    # selected reduced template set, for which experimental annotation exists
    template.file <- fasta.file <- system.file("extdata", "IMGT_data", "templates", 
                                                "Homo_sapiens_IGH_functional_exon.fasta", 
                                                package = "openPrimeR")
    leader.file <- system.file("extdata", "IMGT_data", "templates",
                              "Homo_sapiens_IGH_functional_leader.fasta", 
                              package = "openPrimeR")
    hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")
    # load templates
    template.df <- openPrimeR::read_templates(template.file, hdr.structure, "|", "GROUP", rm.keywords = "partial")
    # assign leader:
    template.df <- assign_binding_regions(template.df, leader.file, NULL)
    # extend forward binding region by one position for overlap cvg starting at V-gene
    template.df <- adjust_binding_regions(template.df, c(-max(template.df$Allowed_End_fw_initial - template.df$Allowed_Start_fw_initial), 0), c(-50,1)) 
    m <- match(colnames(ref.data[[1]]), template.df$ID)
    my.templates <- template.df[m[!is.na(m)],]
    # do the same analysis with the annealing coverage definition!
    tiller.location <- system.file("extdata", "IMGT_data", "primers", "IGHV", 
                            "Tiller2008_1st.fasta", package = "openPrimeR")
    open.location <- system.file("extdata", "IMGT_data", "primers", "IGHV", 
                            "openPrimeR2017.fasta", package = "openPrimeR")
    # create settings for individual tool data generation:
    filename <- system.file("extdata", "settings", 
                "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
    settings <- openPrimeR::read_settings(filename)
    # prevent other binding events:
    conOptions(settings)$allowed_other_binding_ratio <- 0.0
    # allow input nbr of mismatches
    conOptions(settings)$allowed_mismatches <- 12
    # require only intersection with target region
    conOptions(settings)$allowed_region_definition <- "any"
    # compute data
    # activate all 'interesting' cvg constraints as features, but don't filter! (set boundaries high)
    cvg_constraints(settings) <- list("annealing_DeltaG" = c("max" = 0), "primer_efficiency" = c("max" = 1),
                                      "terminal_mismatch_pos" = c("min" = 0))
    tiller.settings <- settings
    PCR(tiller.settings)$annealing_temp <- 57
    tiller.data <- get_tool_primer_data(tiller.location, my.templates, tiller.settings)
    open.settings <- settings
    PCR(open.settings)$annealing_temp <- 55
    open.data <- get_tool_primer_data(open.location, my.templates, open.settings)
    tool.data <- list("Tiller" = tiller.data, "openPrimeR" = open.data)
    feature.data <- vector("list", length(tool.data))
    for (i in seq_along(tool.data)) {
        ident <- names(tool.data)[i]
        primer.df <- tool.data[[i]]$Primers
        template.df <- tool.data[[i]]$Templates
        ref.cvg <- ref.data[[i]]
        data.matrix <- prepare_learning_data(primer.df, template.df, ref.cvg, settings, cvg.type = "basic")
        feature.data[[i]] <- data.matrix
    }
    feature.matrix <- do.call(rbind, feature.data)
    # important: set position of 3' terminus mismatch to primer length if NA (no mismatch)
    pos.mod <- feature.matrix$Position_3terminus 
    # add prior knowledge to positions representation (only consider 3' hexamer)
    pos.mod[is.na(pos.mod) | pos.mod >= 7] <- 7 # no mismatch in 3' hexamer
    pos.mode <- abs(pos.mod - 7)
    # check:
    #feature.matrix$Position_3terminus[pos.mode == 0]
    feature.matrix$Position_3terminusLocal <- pos.mode
    feature.matrix$Position_3terminus[is.na(feature.matrix$Position_3terminus)] <- max(feature.matrix$Position_3terminus, na.rm = TRUE) + 1
    return(feature.matrix)
}
color.to.class <- function(ref.df) {
    cnam <- colnames(ref.df)
    red <- "ff0000"
    orange <- "ed7d31"
    fail.colors  <- c(red, orange)
    conv <- function(color) {
        out <- rep("Covered", length(color))
        idx <- which(color %in% fail.colors) # red-colored -> not amplified
        if (length(idx) != 0) {
            out[idx] <- "Uncovered"
        }
        return(out)
    }
    df <- data.frame(apply(ref.df, 2, conv))
    colnames(df) <- cnam
    return(df)
}

prepare_learning_data <- function(primer.df, template.df, ref.cvg, settings, cvg.type = c("constrained", "basic")) {
    #######
    # prepare tool results for learning:
    #######
    cvg.type <- match.arg(cvg.type)
    full.df <- prepare_mm_plot(primer.df, template.df)
    full.df <- full.df[full.df$Coverage_Type == cvg.type,]
    # select unique binding event for all mismatch contacts:
    df <- plyr::ddply(full.df, c("Primer", "Template", "Group"), plyr::summarize,
                        Number_of_mismatches = unique(substitute(Number_of_mismatches)),
                        primer_efficiency = unique(substitute(primer_efficiency)),
                        annealing_DeltaG = unique(substitute(annealing_DeltaG)),
                        Position_3terminus = min(substitute(Position_3terminus)),
                        All_mismatches = paste(substitute(Position_3prime), collapse = ","))

    # add 3' hexamer mismatch features: 1 if mismatch, 0 if not, per position in the 3' hexamer
    hexa.df <- data.frame(matrix(rep(0, 6 * nrow(df)), nrow = nrow(df), ncol = 6))
    colnames(hexa.df) <- paste0("Mismatch_pos_", seq(6, 1)) # changed labels for considering position starting in the sequence direction 
    #colnames(hexa.df) <- paste0("Mismatch_pos_", seq(1, 6)) # changed labels for considering position starting in the sequence direction 
    mm <- lapply(strsplit(df$All_mismatches, split = ","), function(x) ifelse(x == "NA", NA, as.numeric(x)))
    for (i in seq_along(mm)) {
        cur.mm <- mm[[i]]
        if (!is.na(cur.mm[1])) {
            sel <- which(cur.mm <= 6)
            if (length(sel) != 0) {
                hexa.mm <- cur.mm[sel]
                hexa.df[i, hexa.mm] <- 1
            }
        }
    }
    # take position from the end or just the position in the order of the seq?
    hexa.df <- hexa.df[, rev(seq_len(ncol(hexa.df)))]
    # determine count of hexamer mismatches
    hexa.counts <- apply(hexa.df, 1, sum)
    df <- cbind(df, hexa.df, "Hexamer_Mismatch_Count" = hexa.counts)
    #######
    # prepare experimental data 
    #########
    # turn colors to cvg class:
    ref.cvg <- cbind("Primer" = ref.cvg$Primer, color.to.class(ref.cvg[, 2:(ncol(ref.cvg))]))
    # exclude from ref cvg the templates that are not in template.df!!!!
    my.m <- c(TRUE, sapply(colnames(ref.cvg)[-1], function(x) length(grep(x, template.df$ID, fixed = TRUE)) != 0))
    ref.cvg <- ref.cvg[, my.m]
    # TODO:: treatment of replicates -> include all / keep consensus / something else? for now: only assign consensus
    ref.cvg.consensus <- plyr::ddply(ref.cvg, "Primer", plyr::catcolwise(function(x) ifelse(all(x == "Covered"), "Covered", "Uncovered")))
    # convert from wide to long format
    ref.cvg.matrix <- reshape2::melt(ref.cvg.consensus, "Primer", value.name = "Experimental_Coverage", variable.name = "Template")
    # assign labels to the learning data
    combi.df <- merge(df, ref.cvg.matrix, by=c("Primer","Template")) # NA's match
    if (nrow(combi.df) != nrow(df)) {
        warning("combi.df excluded some obsevations; probably because number of allowed mismatches was too low when evaluating the primer coverage -> set it higher for creating a full training set.")
        print(paste("Obs. in tool data frame: ", nrow(df)))
        print(paste("Obs. in merged data frame: ", nrow(combi.df)))
    }
    # assign features from evaluation of primers:
    sel.features <- names(constraints(settings))
    excl.features <- c("cross_dimerization", "melting_temp_diff") # not single-primer based constraints
    sel.features <- setdiff(sel.features, excl.features)
    # rename some features for matching:
    sel.features[sel.features == "melting_temp_range"] <- "melting_temp"
    sel.features[sel.features == "secondary_structure"] <- "Structure_deltaG"
    sel.features[sel.features == "self_dimerization"] <- "Self_Dimer_DeltaG"
    feature.idx <- sapply(sel.features, function(x) grep(paste0("^",x), colnames(primer.df))[1])
    feature.df <- cbind(Primer = primer.df$ID, primer.df[, feature.idx])
    # extract features and assign to the combined data frame
    augmented.df <- merge(combi.df, feature.df, by = "Primer")
    # store info about primer set in rownames
    rownames(augmented.df) <- paste0(rownames(augmented.df), "_", unique(primer.df$Run))
    augmented.df$Experimental_Coverage <- factor(augmented.df$Experimental_Coverage, levels = c("Uncovered", "Covered"))
    # add run
    augmented.df$Run <- unique(primer.df$Run)
    return(augmented.df)
}
devtools::load_all("src/openPrimeR")
######
# get reference coverage data frame:
######
xls.file <- "data/PCR_ref_data.xlsx"
ref.data <- get_ref_data (xls.file)
######
# get feature matrix for supervised learning
#####
feature.matrix <- get_learning_matrix(ref.data)
#########
# Store the reference data, feature matrix:
#########
out.loc <- file.path(system.file("data",  package = "openPrimeR"), "RefCoverage.rda")
save(ref.data, feature.matrix, file = out.loc, compress = "xz")

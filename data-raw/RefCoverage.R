#########
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
get_ref_data <- function(cvg.folders) {
    #############################
    consensus.fun <- function(x) {
        count1 <- length(which(x == 1))
        count0 <- length(which(x == 0))
        if (count1 > count0) {
            return(1)
        } else {
            return(0)
        }
    }
    agreement.rate <- function(x) {
        count1 <- length(which(x == 1))
        count0 <- length(which(x == 0))
        majority <- ifelse(count1 > count0, 1, 0)
        agreement.rate <- length(which(x == majority)) /length(x)
        return(agreement.rate)
    }
    library(plyr)
    library(reshape2)
    library(ggplot2)
    ignore.cols <- c("Evaluator", "experiment")
    ref.calls <- vector("list", length(cvg.folders))
    for (i in seq_along(cvg.folders)) { # primer sets
        csvs <- list.files(cvg.folders[i], full.names = TRUE)
        set.name <- basename(cvg.folders[i])
        ref.results <- vector("list", length(csvs))
        for (j in seq_along(csvs)) { # reviewer results
            csv <- read.csv(csvs[j], check.names = FALSE)
            csv.j <- csv[, !colnames(csv) %in% ignore.cols]
            # determine consensus across individual experiments
            res <- ddply(csv.j, c("Primer"), numcolwise(consensus.fun))
            res <- cbind("Evaluator" = unique(csv$Evaluator), res)
            ref.results[[j]] <- res
        }
        cvg.matrix <- do.call(rbind, ref.results)
        cvg.matrix.j <- cvg.matrix[, !colnames(cvg.matrix) %in% ignore.cols]
        # compute agreement rate per primer template for all the reviewers from this
        agree.df <- ddply(cvg.matrix.j, c("Primer"), numcolwise(agreement.rate))
        plot.df <- melt(agree.df, variable.name = "Template", value.name = "Agreement")
        # plot 
        p <- ggplot() + geom_boxplot(data = plot.df, aes(x = Primer, y = Agreement))
        ggsave(paste0("agreement_", set.name, ".png"), p)
        ref.summary <- ddply(cvg.matrix.j, c("Primer"), numcolwise(consensus.fun))
        ref.calls[[i]] <- ref.summary
    }
    # gather results for Tiller and openPrimeR in a data frame
    result.df <- do.call(rbind, ref.calls)
    # convert to machine-readable format
    res <- melt(result.df, variable.name = "Template", value.name = "Experimental_Coverage")
    res$Experimental_Coverage <- factor(ifelse(res$Experimental_Coverage == 1, "Amplified", "Unamplified"), levels = c("Unamplified", "Amplified"))
    return(res)
}
get_tool_primer_data <- function(primer.location, template.df, cur.settings) {
    # load primers
    primer.df <- openPrimeR::read_primers(primer.location, "_fw", "_rev")
    primer.df$ID <- factor(gsub("_fw", "", as.character(primer.df$ID)))
    primer.df <- openPrimeR::check_constraints(primer.df, template.df, cur.settings)
    template.df <- openPrimeR::update_template_cvg(template.df, primer.df)
    out <- list("Primers" = primer.df, "Templates" = template.df)
    return(out)
}
get_ref_data_old <- function(xls.file) {
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
get_learning_matrix <- function(ref.data, template.df, tiller.location, open.location) {
    # primer set to use for calibration:
    mode.directionality <- "fw"
    # create tool data for Tiller / openPrimeR:
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
    tiller.data <- get_tool_primer_data(tiller.location, template.df, tiller.settings)
    open.settings <- settings
    PCR(open.settings)$annealing_temp <- 55
    open.data <- get_tool_primer_data(open.location, template.df, open.settings)
    tool.data <- list("Tiller" = tiller.data, "openPrimeR" = open.data)
    feature.data <- vector("list", length(tool.data))
    for (i in seq_along(tool.data)) {
        ident <- names(tool.data)[i]
        primer.df <- tool.data[[i]]$Primers
        template.df <- tool.data[[i]]$Templates
        data.matrix <- prepare_learning_data(primer.df, template.df, ref.data, settings, cvg.type = "basic")
        feature.data[[i]] <- data.matrix
    }
    feature.matrix <- do.call(rbind, feature.data)
    # encode mismatches locally: hexamer encoding
    pos.mod <- feature.matrix$Position_3terminus 
    pos.mod[is.na(pos.mod) | pos.mod >= 7] <- 7 # no mismatch in 3' hexamer
    pos.mod <- abs(pos.mod - 7)
    feature.matrix$Position_3terminusLocal <- pos.mod
    #feature.matrix$Position_3terminus[is.na(feature.matrix$Position_3terminus)] <- max(feature.matrix$Position_3terminus, na.rm = TRUE) + 1
    feature.matrix$Number_of_mismatches_hexamer <- feature.matrix$Mismatch_pos_1 +  feature.matrix$Mismatch_pos_2 +  
                                                    feature.matrix$Mismatch_pos_3 +  feature.matrix$Mismatch_pos_4 +  
                                                    feature.matrix$Mismatch_pos_5 +  feature.matrix$Mismatch_pos_6
    return(feature.matrix)
}

get_learning_matrix_old <- function(ref.data) {
    data(Ippolito) # load templates
    # primer set to use for calibration:
    mode.directionality <- "fw"
    # create tool data for Tiller / openPrimeR:
    # selected reduced template set, for which experimental annotation exists
    template.file <- fasta.file <- system.file("extdata", "IMGT_data", "templates", 
                                                "Homo_sapiens_IGH_functional_exon_exp.fasta", 
                                                package = "openPrimeR")
    leader.file <- system.file("extdata", "IMGT_data", "templates",
                              "Homo_sapiens_IGH_functional_leader_exp.fasta", 
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
    feature.matrix$Position_3terminusLocal <- pos.mode
    feature.matrix$Position_3terminus[is.na(feature.matrix$Position_3terminus)] <- max(feature.matrix$Position_3terminus, na.rm = TRUE) + 1
    feature.matrix$Number_of_mismatches_hexamer <- feature.matrix$Mismatch_pos_1 +  feature.matrix$Mismatch_pos_2 +  
                                                    feature.matrix$Mismatch_pos_3 +  feature.matrix$Mismatch_pos_4 +  
                                                    feature.matrix$Mismatch_pos_5 +  feature.matrix$Mismatch_pos_6
    return(feature.matrix)
}
color.to.class <- function(ref.df) {
    cnam <- colnames(ref.df)
    red <- "ff0000"
    orange <- "ed7d31"
    fail.colors  <- c(red, orange)
    conv <- function(color) {
        out <- rep("Amplified", length(color))
        idx <- which(color %in% fail.colors) # red-colored -> not amplified
        if (length(idx) != 0) {
            out[idx] <- "Unamplified"
        }
        return(out)
    }
    df <- data.frame(apply(ref.df, 2, conv))
    colnames(df) <- cnam
    return(df)
}

prepare_learning_data <- function(primer.df, template.df, ref.data, settings, cvg.type = c("constrained", "basic")) {
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
                        All_mismatches = paste(substitute(Position_3prime), collapse = ","),# mismatches relative to the 3' end (pos 1 indicates 3' terminus)
                        Worst_Case_MismatchLocal = paste(substitute(Position_3terminusLocal), collapse = ","))  # worst case mismatches where 1 -> first pos in 3' hexamer, 6 -> 3' terminus
                        #Binding_Pos_Start = unique(substitute(Binding_Position_Start_fw)))
    # add some primer features
    m <- match(df$Primer, primer.df$ID)
    df$Primer_Sequence <- primer.df$Forward[m]
    df$Terminal_Dinucleotide <- substr(primer.df$Forward[m], nchar(primer.df$Forward[m]) - 1, nchar(primer.df$Forward[m]))
    # add some template features; TODO: need to adjust prepare_mm_plot to output the binding stretch or binding position in the templates
    m <- match(df$Template, template.df$ID)
    template.seqs <- template.df$Sequence[m]
    # add 3' hexamer mismatch features: 1 if mismatch, 0 if not, per position in the 3' hexamer: 
    # pos 1 -> first position in 3' hexamer, pos6: last position in 3' hexamer (<=> 3' terminus)
    hexa.df <- data.frame(matrix(rep(0, 6 * nrow(df)), nrow = nrow(df), ncol = 6))
    colnames(hexa.df) <- paste0("Mismatch_pos_", seq(1, 6)) # changed labels for considering position starting in the sequence direction 
    mm <- lapply(strsplit(df$All_mismatches, split = ","), function(x) ifelse(x == "NA", NA, as.numeric(x)))
    for (i in seq_along(mm)) {
        cur.mm <- mm[[i]]
        if (!is.na(cur.mm[1])) {
            sel <- which(cur.mm <= 6)
            if (length(sel) != 0) {
                hexa.mm <- abs(cur.mm[sel] - 7)
                hexa.df[i, hexa.mm] <- 1
            }
        }
    }
    # determine count of hexamer mismatches
    hexa.counts <- apply(hexa.df, 1, sum)
    df <- cbind(df, hexa.df, "Hexamer_Mismatch_Count" = hexa.counts)
    #######
    # prepare experimental data 
    #########
    # assign labels to the learning data
    combi.df <- merge(df, ref.data, by=c("Primer","Template")) # NA's match
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
    augmented.df$Experimental_Coverage <- factor(augmented.df$Experimental_Coverage, levels = c("Unamplified", "Amplified"))
    # add run
    augmented.df$Run <- unique(primer.df$Run)
    return(augmented.df)
}
devtools::load_all("src/openPrimeR") # load openPrimeR package
# get reference coverage for primer sets
cvg.folder <- system.file("data-raw", "coverage_data", package = "openPrimeR")
cvg.folders <- c(file.path(cvg.folder, "Tiller"), file.path(cvg.folder, "openPrimeR"))
ref.data <- get_ref_data(cvg.folders)
#save(ref.data, file = out.loc, compress = "xz")
# load templates
template.file <- file.path(cvg.folder, "IGHV_templates.fasta") # NB: use these templates to assign coverage!!! (TODO: check leader ... other file from christoph?)
# define 5' UTR + leader as allowed binding region
allowed.file <- file.path(cvg.folder, "IGHV_allowed_region.fasta")
template.df <- read_templates(template.file, c("GROUP"), delim = "|") 
template.df <- assign_binding_regions(template.df, fw = allowed.file, rev = NULL)
# extend forward binding region by one position in order to set coverage definition where any overlap with binding region is tolerated (ensure that primers bind at the start position at the latest)
template.df <- adjust_binding_regions(template.df, c(-max(template.df$Allowed_End_fw_initial - template.df$Allowed_Start_fw_initial), 0), c(-50,1)) 
# load primers from adjusted paths due to adapted filenames fitting for reference primer IDs
tiller.location <- file.path(cvg.folder, "Tiller2008_1st.fasta")
open.location <- file.path(cvg.folder,"openPrimeR2017.fasta")
######
# get feature matrix for supervised learning
#####
feature.matrix <- get_learning_matrix(ref.data, template.df, tiller.location, open.location)
out.loc <- file.path(system.file("data",  package = "openPrimeR"), "RefCoverage.rda")
#########
# Store the reference data, feature matrix:
#########
save(ref.data, feature.matrix, file = out.loc, compress = "xz")
############
if (FALSE) {
    # old code:
    xls.file <- system.file("data-raw", "PCR_ref_data.xlsx", package = "openPrimeR")
    hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")
    template.df <- read_templates(system.file("extdata", "IMGT_data", "templates", 
                                                    "Homo_sapiens_IGH_functional_exon_exp.fasta", 
                                                    package = "openPrimeR"), hdr.structure, "|")
    template.df <- assign_binding_regions(template.df, fw = system.file("inst", "extdata", "IMGT_data", "templates", "Homo_sapiens_IGH_functional_leader_exp.fasta", package = "openPrimeR"), rev = NULL)
    ref.data <- get_ref_data (xls.file)
    feature.matrix <- get_learning_matrix(ref.data)
}


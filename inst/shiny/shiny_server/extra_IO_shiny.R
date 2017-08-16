################
# Input/Output functions for Shiny App
#################

view.input.sequences <- function(template.df) {
    # shiny output after reading template sequences
    if (nrow(template.df) == 0 || length(template.df) == 0) {
        return(NULL)
    }
    template.df <- asS3(template.df) # modifying columns -> type can't be preserved
    excl.col <- c("Header", "Group", "Identifier", "Sequence_Length", "Allowed_Start_fw", 
        "Allowed_End_fw", "Allowed_Start_rev", "Allowed_End_rev", 
        "Run", "Allowed_Start_fw_initial", "Allowed_End_fw_initial", 
        "Allowed_Start_rev_initial", "Allowed_End_rev_initial", "InputSequence")
    # add all "ali" columns to exclusion:
    ali.cols <- colnames(template.df)[grep("_ali", colnames(template.df))]
    excl.col <- c(excl.col, ali.cols)
    view.df <- openPrimeR:::exclude.cols(excl.col, template.df)
    # remove all columns where all values are missing
    excl.idx <- which(unlist(lapply(seq_len(ncol(view.df)), function(x) all(view.df[,x] == "" | is.na(view.df[,x]) | view.df[,x] == 
        "NA-NA"))))
    if (length(excl.idx) != 0) {
        view.df <- view.df[, -excl.idx]
    }
    view.df <- openPrimeR:::modify.col.rep(view.df)
    return(view.df)
}

update.sample.name <- function(df, sample.name) {
    # updates the Run column of a data frame
    if (length(df) == 0 || nrow(df) == 0) {
        return(NULL)
    }
    if (length(sample.name) == 0) {
        return(df)
    } else {
        df$Run <- sample.name
    }
    return(df)
}
view.lex.sequences <- function(template.df) {
    # shiny output after leaders have been assigned
    if (nrow(template.df) == 0 || length(template.df) == 0) {
        return(NULL)
    }
    template.df <- asS3(template.df)
    excl.col <- c() # nothing to exclude here at the moment
    view.df <- openPrimeR:::exclude.cols(excl.col, template.df)
    # show value ranges in dash notation as one column only show if we have some
    # annotated leaders..
    if ("Allowed_Start_fw" %in% colnames(view.df)) {
        if (any(!is.na(template.df$Allowed_Start_fw))) {
            leader.range.fw <- paste0(template.df$Allowed_Start_fw, "-", template.df$Allowed_End_fw)
            view.df$Allowed_Start_fw <- leader.range.fw
        }
    }
    if ("Allowed_Start_rev" %in% colnames(view.df)) {
        if (any(!is.na(template.df$Allowed_Start_rev))) {
            leader.range.rev <- paste0(template.df$Allowed_Start_rev, "-", template.df$Allowed_End_rev)
            view.df$Allowed_Start_rev <- leader.range.rev
        }
    }
    change.cols <- c("Allowed_Start_fw", "Allowed_Start_rev")
    new.names <- c("Allowed Binding Range (fw)", "Allowed Binding Range (rev)")
    for (i in seq_along(change.cols)) {
        colnames(view.df)[colnames(view.df) == change.cols[i]] <- new.names[i]
    }
    view.df <- view.input.sequences(view.df)
    return(view.df)
}

view.mismatch.table <- function(mismatch.table) {
    if (length(mismatch.table) == 0 || nrow(mismatch.table) == 0) {
        return(NULL)
    }
    excl.col <- c("NT_nbr_mm", "AA_nbr_mm", "Seq_NT", "Primer_NT", "Seq_AA", "Primer_AA", 
        "Primer_ID", "Primer_Seq")
    view.df <- openPrimeR:::exclude.cols(excl.col, mismatch.table)
    col.order <- c("Primer", "Template", "Alignment_NT", "Mutation_Type", 
        "Alignment_AA", "Comment")
    view.df <- view.df[, col.order]
    view.df <- openPrimeR:::modify.col.rep(view.df)
    return(view.df)
}
view.cvg.sequences <- function(template.df, primer.df) {
    # viewing templates with annotated coverage
    if (nrow(template.df) == 0 || length(template.df) == 0) {
        return(NULL)
    }
    if (!"primer_coverage" %in% colnames(template.df)) {
        # nothing to change here
        return(view.lex.sequences(template.df))
    }
    template.df <- asS3(template.df)
    excl.col <- c() # no columns to exclude
    view.df <- openPrimeR:::exclude.cols(excl.col, template.df)
    # convert from identifier coverage to ID coverage
    cvg.cols <- c("Covered_By_Primers", "Covered_By_Primers_fw", 
                    "Covered_By_Primers_rev")
    for (i in seq_along(cvg.cols)) {
        col <- cvg.cols[i]
        if (col %in% colnames(view.df)) {
            view.df[, col] <- unlist(openPrimeR:::covered.primers.to.ID.string(view.df[, col], primer.df)) 
       }
    }
    # re-order
    col.order <- c("ID", "Group", 
        "primer_coverage", "primer_coverage_fw", 
        "primer_coverage_rev", "Allowed_fw", "Allowed_rev",
        "Allowed_Start_fw", "Allowed_End_fw",
        "Allowed_Start_rev", "Allowed_End_rev",
        "Covered_By_Primers", "Covered_By_Primers_fw",
        "Covered_By_Primers_rev"
        )
    other.cols <- setdiff(colnames(view.df), col.order)  # don't care about these cols in ordering
    view.df <- view.df[, c(col.order, other.cols)]
    # don't order by cvg -> misleading
    #view.df <- view.df[order(view.df$primer_coverage, decreasing = TRUE), ]
    # remove columns of unrelevant direction
    mode.directionality <- openPrimeR:::get.analysis.mode(primer.df)
    # special.cols: redundant for single direction
    special.cols <- c("Covered_By_Primers", "primer_coverage")
    excl <- NULL
    if (mode.directionality == "fw") {
        excl <- grep("_rev", colnames(view.df))
        idx <- which(colnames(view.df) %in% paste0(special.cols, "_fw"))
        excl <- c(excl, idx)
    } else if (mode.directionality == "rev") {
        excl <- grep("_fw", colnames(view.df))
        idx <- which(colnames(view.df) %in% paste0(special.cols, "_rev"))
        excl <- c(excl, idx)
    }
    if (length(excl) != 0) {
        view.df <- view.df[, -excl]
    }
    view.df <- view.lex.sequences(view.df)
    return(view.df)
}
view.subset.primers <- function(primer.df, template.df, mode.directionality, view.cvg.individual = "inactive") {
    if (length(primer.df) == 0) {
        return(NULL)
    } else if (nrow(primer.df) == 0) {
        return(primer.df)
    }
    view.df <- asS3(primer.df)
    # convert covered template seqs to group representation
    cvd <- openPrimeR:::covered.seqs.to.ID.string(as.character(view.df$Covered_Seqs), template.df)
    if (length(unique(template.df$Group)) >= 2 && view.cvg.individual == "inactive") {
        # show gene group instead of identifiers
        idx <- openPrimeR:::covered.seqs.to.idx(as.character(view.df$Covered_Seqs), template.df)
        cvd <- openPrimeR:::string.list.format(sapply(seq_along(idx), function(x) paste(template.df[idx[[x]], 
            "Group"], collapse = ",")))
    }
    view.df$Covered_Seqs <- unlist(cvd)
    # percent format the cvg ratio
    if ("Coverage_Ratio" %in% colnames(view.df)) {
        view.df$Coverage_Ratio <- paste(round(view.df$Coverage_Ratio, 4) * 100, "%", 
            sep = "")
    }
    excl <- NULL
    if (mode.directionality == "fw") {
        excl <- grep("_rev", colnames(view.df))
    } else if (mode.directionality == "rev") {
        excl <- grep("_fw", colnames(view.df))
    }
    excl.cols <- c("Direction", "primer_length_fw", "primer_length_rev")
    excl.idx <- sapply(excl.cols, function(x) grep(x, colnames(view.df)))
    excl.idx <- unique(c(excl.idx, excl))
    if (length(excl.idx) != 0) {
        view.df <- view.df[, -excl.idx]
    }
    view.df <- openPrimeR:::view.input.primers(view.df, mode.directionality)
    return(view.df)
}
view.filtered.primers.all <- function(primer.df, template.df, mode.directionality, view.cvg.individual) {
    # view evaluated primers (name is misleading)
    if (length(primer.df) == 0) {
        return(NULL)
    } else if (nrow(primer.df) == 0) {
        return(primer.df)
    }
    primer.df <- asS3(primer.df)
    if ("primer_coverage" %in% colnames(primer.df)) {
        view.df <- openPrimeR:::view.cvg.primers(primer.df, template.df, mode.directionality, view.cvg.individual)
    } else {
        view.df <- primer.df
    }
    return(view.df)
}
view.filtered.primers <- function(primer.df, template.df, mode.directionality, view.cvg.individual) {
    # view excluded primers message('viewing filtered primers:') message(primer.df)
    if (length(primer.df) == 0) {
        return(NULL)
    } else if (nrow(primer.df) == 0) {
        # still show empty df
        return(primer.df)
    }
    primer.df <- asS3(primer.df)
    excl.col <- "constraints_passed"
    # also remove all of the eval columns:
    eval.cols <- colnames(primer.df)[grep("EVAL_", colnames(primer.df))]
    excl.col <- c(excl.col, eval.cols)
    view.df <- openPrimeR:::exclude.cols(excl.col, primer.df)
    view.df <- view.filtered.primers.all(view.df, template.df, mode.directionality, view.cvg.individual)
    return(view.df)
}
view.optimized.primers <- function(primer.df, template.df, mode.directionality, view.cvg.individual) {
    if (length(primer.df) == 0) {
        return(NULL)
    } else if (nrow(primer.df) == 0) {
        # still display an empty df
        return(primer.df)
    }
    view.df <- asS3(primer.df)
    if ("Cumulative_Coverage_Ratio" %in% colnames(view.df)) {
        view.df$Cumulative_Coverage_Ratio <- paste(round(view.df$Cumulative_Coverage_Ratio, 
            4) * 100, "%", sep = "")
    }
    view.df <- view.filtered.primers(view.df, template.df, mode.directionality, view.cvg.individual)
    # reorder columns: put ID first
    c.order <- c("ID")
    c.o <- setdiff(colnames(view.df), c.order)  # don't care about these cols in ordering
    cols <- c(c.order, c.o)
    view.df <- view.df[, cols]
    return(view.df)
}
write.out.constraints <- function(constraint.settings, file.out) {
    # TODO: Deprecated?
    write("Constraint settings overview", file.out)  # overwrite file on new call
    for (i in seq_along(constraint.settings)) {
        name <- names(constraint.settings)[i]
        val <- constraint.settings[[i]]
        val.names <- names(val)
        val.string <- paste(val.names, "=", val)
        text <- paste(name, ": ", paste(val.string, collapse = ", ", sep = ""), sep = "")
        write(text, file.out, append = TRUE, ncolumns = 1000)
    }
}

constraint.file.to.list <- function(constraint.files) {
    constraints <- vector("list", length(constraint.files))
    for (i in seq_along(constraint.files)) {
        if (is.na(constraint.files[i])) {
            filtering.constraints <- list()
        } else {
            data <- read_settings(constraint.files[i])
            # only the actually applied constraints interest us
            filtering.constraints <- openPrimeR:::filtersUsed(data)
        }
        constraints[[i]] <- filtering.constraints
    }
    return(constraints)
}

comparison.constraint.equivalences <- function(constraints) {
    M <- matrix(rep(FALSE, length(constraints) * length(constraints)), nrow = length(constraints))
    for (i in seq_along(constraints)) {
        for (j in seq_along(constraints)) {
            # message(paste(i, ',', j, sep = ''))
            M[i, j] <- compare.constraints(constraints[[i]], constraints[[j]])
        }
    }
    # iterate through M to find equivalence classes
    result <- vector("list", length(constraints))  # worst-case length
    nbr.classes <- 0
    sets.covered <- vector("list", length(constraints))
    for (i in seq_along(constraints)) {
        idx <- which(M[i, ])  # columns (primer sets) that have the same constraints
        if (!all(idx %in% unlist(sets.covered))) {
            # sets.covered <- c(sets.covered, idx)
            if (length(idx) != 0) {
                nbr.classes <- nbr.classes + 1
                sets.covered[[nbr.classes]] <- idx
                result[[nbr.classes]] <- constraints[[i]]
            }
        }
    }
    result <- result[1:nbr.classes]
    sets.covered <- sets.covered[1:nbr.classes]
    # order equivalence classes by size
    # such that 1st class is the largest
    l <- sapply(sets.covered, length)
    o <- order(l, decreasing = TRUE)
    sets.covered <- sets.covered[o]
    result <- result[o]
    equivalence.idx <- rep(0, length(constraints))
    for (i in seq_along(sets.covered)) {
        equivalence.idx[sets.covered[[i]]] <- i
    }
    if (length(result) == 1 && is.null(result[[1]])) {
        # no constraints available
        out <- list(Index = 1, Constraints = result)
    } else {
        out <- list(Index = equivalence.idx, Constraints = result)
    }
    return(out)
}
change.extension <- function(x, ext) {
    out <- sub("\\.[[:alnum:]]+$", "", x)
    out <- paste(out, ".", ext, sep = "")
    return(out)
}
call.IMGT.settings.script <- function(imgt.settings.location) {
    # calls the IMGT settings script to generate output files for 'get.IMGT.settings'
    script.location <- file.path(server.src.folder, "extra_IMGT_template_options.py")
    # args: phantomjs, out-folder
    phantomjs.loc <- Sys.which("phantomjs")
    if (phantomjs.loc == "") {
        stop("PhantomJS not available.")
    }
    args <- paste(normalizePath(phantomjs.loc), imgt.settings.location, sep = " ")
    call <- paste(script.location, " ", args, sep = "")
    #message(call)
    # store all IMGT field options in text files
    system(call, ignore.stdout = TRUE)  
}
get.IMGT.settings <- function() {
    # use a python script to retrieve possible IMGT input options for shiny frontend
    # options are stored in files and are read into R
    imgt.settings.location <- file.path(app.data.folder, "IMGT_options")
    if (!dir.exists(imgt.settings.location)) {
        # only retrieve options from IMGT when they aren't available yet
        dir.create(imgt.settings.location)
        call.IMGT.settings.script(imgt.settings.location)
    }
    # retrieve data from files
    files <- list.files(imgt.settings.location)
    options <- vector("list", length(files))
    names(options) <- files
    for (i in seq_along(files)) {
        f <- file.path(imgt.settings.location, files[i])
        data <- readLines(f)
        options[[i]] <- data
    }
    return(options)
}
get.supplied.comparison.template.path <- function(locus) {
    # loads locally stored csv containing template analysis results
    if (length(locus) == 0 || locus == "") {
        return(NULL)
    }
    fname <- paste(locus, "_templates.csv", sep = "")
    load.path <- system.file("extdata", "IMGT_data", "comparison", "templates", 
                 fname, package = "openPrimeR")
    if (!file.exists(load.path)) {
        warning(paste("Template comparison data for ", locus, " not found! Specified comparison path was: ", 
            load.path, sep = ""))
        return(NULL)
    }
    res <- list(datapath = load.path, name = locus)
    return(res)
}
primer.set.choices <- function(primers) {
    if (length(primers) == 0) {
        return(NULL)
    }
    paths <- primers
    locus <- strsplit(primers[1], split = "/")[[1]]
    locus <- locus[length(locus) - 1]
    # get basenames
    primers <- basename(primers)
    # remove 'IPS' annotation
    idx <- grep("IPS", primers)
    mod <- unlist(lapply(strsplit(primers[idx], split = "_"), function(x) x[[2]]))
    primers[idx] <- mod
    # remove extension
    primers <- sub("^([^.]*).*", "\\1", primers)
    # remove name of locus from set name
    idx <- grep(locus, primers)
    mod <- sub(paste("_", locus, sep = ""), "", primers[idx])
    primers[idx] <- mod
    o <- order(primers)  # order by author
    primers <- primers[o]
    paths <- paths[o]  # order paths also ..
    names(paths) <- primers
    return(paths)
}
comparison.primer.choices <- function(locus) {
    if (length(locus) == 0) {
        return(NULL)
    }
    path <- system.file("extdata", "IMGT_data", "comparison",
                        "primer_sets", package = "openPrimeR")
    path <- file.path(path, locus)
    if (!dir.exists(path)) {
        warning(paste("Primer comparison folder ", path, " not found!", sep = ""))
        return(NULL)
    }
    fnames <- list.files(path, pattern = "*.csv", full.names = TRUE)
    primers <- basename(fnames)
    # remove 'IPS' annotation
    idx <- grep("IPS", primers)
    mod <- unlist(lapply(strsplit(primers[idx], split = "_"), function(x) x[[3]]))
    primers[idx] <- mod
    # remove extension
    primers <- gsub("^([^.]*).*", "\\1", primers)
    # remove name of locus from set name
    idx <- grep(locus, primers)
    mod <- sapply(strsplit(primers[idx], split = "_"), function(x) paste(x[!grepl(locus, 
        x)], collapse = "_"))
    primers[idx] <- mod
    o <- order(primers)  # order by author
    out <- fnames[o]
    names(out) <- primers[o]
    return(out)
}
retrieve.IMGT.templates <- function(species, locus, func, refresh, rm.partials = FALSE) {
    # species, locus, function: parameters for IMGT data retrieval refresh:
    # TRUE/FALSE -> should data be retrieved (TRUE) or should stored data be used if
    # available (FALSE)
    # data sets for which we have novel sequences not in IMGT
    new.data.sets <- c("IGH")
    new.data.available <- species == "Homo sapiens" && locus %in% new.data.sets && func == "functional"
    if (new.data.available) {
        # there's new data available -> load these data from package data
        fname.leader <- get.IMGT.fname(species, locus, func, "leader_exp")
        fname.exon <- get.IMGT.fname(species, locus, func, "exon_exp")
        ret <- c(fname.exon, fname.leader)
        if (!all(file.exists(ret))) {
            warning("File did not exist: ", paste0(ret, collapse = ","))
        }
        return(ret)
    } else {
        fname.leader <- get.IMGT.fname(species, locus, func, "leader")
        fname.exon <- get.IMGT.fname(species, locus, func, "exon")
    }
    ret <- c(fname.exon, fname.leader)
    existing.files <- file.exists(ret)
    if (all(existing.files) && !refresh) {
        # don't recompute if results are already available.
        return(ret)
    }
    imgt.script <- file.path(server.src.folder, "extra_IMGT_template_set_extractor.py")
    if (Sys.which("phantomjs") == "") {
        warning("PhantomJS is not in your path. Please install it.")
        return(NULL)
    }
    if (!openPrimeR:::selenium.installed()) {
        warning("Selenium for Python is not available. Please install it.")
        return(NULL)
    }
    base.opts <- paste(normalizePath(Sys.which("phantomjs")), fname.leader, 
        fname.exon, sep = " ")
    input.opts <- paste(paste("\"", species, "\"", sep = ""), paste("\"", locus, 
        "\"", sep = ""), paste("\"", func, "\"", sep = ""), sep = " ")
    cmd <- paste(imgt.script, base.opts, input.opts, sep = " ")
    #print(cmd)
    # warning(cmd) # printout for shiny-server (docker)
    status <- system(cmd, ignore.stdout = TRUE)
    if (status != 0) {
        warning("Failed to retrieve templates from IMGT for unknown reason.")
        return(NULL)
    } else {
        return(ret)
    }
}
get.IMGT.fname <- function(species, locus, func, type) {
    # species, locus, func: args for IMGT db retrieval type: leader or exon replace
    # gaps with underscores for fnames
    species <- gsub(" ", "_", species)
    locus <- gsub(" ", "_", locus)
    func <- gsub(" ", "_", func)
    types <- c("leader", "exon")
    template.folder <- system.file("extdata", "IMGT_data", "templates",
                        package = "openPrimeR")
    names <- file.path(template.folder, paste(species, "_", locus, "_", func, 
        "_", type, ".fasta", sep = ""))
    return(names)
}

get.available.settings <- function(app.settings.folder, taq.PCR = NULL, analysis.mode = NULL, initial = FALSE) {
    # select the setting files that correspond to the user input
    setting.files <- list.files(app.settings.folder, ".xml", full.names = TRUE)
    # select only setting corresponding to selected polymerase
    if (length(taq.PCR) != 0) {
        if (taq.PCR) {
            # taq PCR
            sel <- grep("_Taq_", basename(setting.files))
        } else {
            # non-taq PCR
            sel <- grep("_Non-Taq_", basename(setting.files))
        }
        setting.files <- setting.files[sel]
    }
    if (length(analysis.mode) != 0) {
        if (analysis.mode == "evaluate" || analysis.mode == "compare") {
            sel <- grep("evaluate", basename(setting.files))
        } else { # design
            sel <- grep("design", basename(setting.files))
        }
        setting.files <- setting.files[sel]
    }
    if (length(setting.files) > 1 & initial) {
        # don't choose a setting yet (UI will still change)?
        setting.files <- setting.files[grep("evaluate", setting.files)[1]]
    }
    return(setting.files)
}
get.available.settings.view <- function(app.settings.folder, taq.PCR = NULL, analysis.mode = NULL, initial = FALSE) {
    fnames <- get.available.settings(app.settings.folder, taq.PCR, analysis.mode, initial = initial)
    if (length(fnames) == 0) {
        return(NULL)
    }
    bnames <- basename(fnames)
    bnames <- sub("^([^.]*).*", "\\1", bnames)
    return(bnames)
}

dimer.text.info <- function(dimer.data, primer.df, deltaG.cutoff) {
    # dimer.data: worst-case conformation for primer pairs
    if (length(dimer.data) == 0) {
        return(NULL)
    }
    dimer.idx <- which(dimer.data$DeltaG < deltaG.cutoff)
    possible.dimers <- nrow(dimer.data)  # the number of primer pairs above the score cutoff
    N.dimers <- length(dimer.idx)
    ratio <- round((N.dimers/possible.dimers) * 100, 2)
    ID.col <- ifelse("Primer" %in% colnames(dimer.data), "Primer", c("Primer_1", "Primer_2"))
    N.primers <- length(unique(unlist(dimer.data[dimer.idx, ID.col])))
    ratio.primers <- round(N.primers/nrow(primer.df) * 100, 2)
    text <- paste("At the current &Delta;G cutoff of ", deltaG.cutoff, " kcal/mol, ", 
        N.dimers, " (", ratio, "%) of ", possible.dimers, " primer pairings are considered dimerizing and ", 
        N.primers, " (", ratio.primers, "%) of ", nrow(primer.df), " primers are involved in these dimerizations.", 
        sep = "")
    return(text)
}

create.design.string <- function(allowed.mismatches, mode.directionality, 
                                 init.mode, opti.algo, template.df,
                                 required.cvg) {
    init <- ifelse(init.mode == "naive", "divergent", "related")
    algo <- ifelse(opti.algo == "Greedy", "a greedy algorithm",
                    "an integer linear program")
    mm.string <- ifelse(allowed.mismatches <= 1, "mismatch", "mismatches")
    dir <- ifelse(mode.directionality == "fw", "forward",
                ifelse(mode.directionality == "rev", "reverse",
                        "pairs of"))
    req.cvg <- paste0(round(required.cvg * 100, 2), "%")
    msg <- paste0("Design ", dir, " primers for ", nrow(template.df), " ", init, " template sequences with at most ", allowed.mismatches, " ", mm.string, " using ", algo, " (target coverage: ", req.cvg, ")?")
    return(msg)
}

myHeaderPanel <- function (title, windowTitle = title, style = NULL) 
{
    tagList(tags$head(tags$title(windowTitle)), 
        div(style = NULL, id = "headerPanel", class = "col-sm-12", title)
    )
}


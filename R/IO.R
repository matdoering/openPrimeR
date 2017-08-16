################
# Input/Output
################

#' Coverage Info Text
#'
#' Creates a string with information on the coverage
#'
#' @param stats Data frame with coverage statistics.
#' @param selected.group Retrieve information for a subgroup of templates only.
#' @param ident An identifier for the coverage.
#' @return A string with information on the coverage.
#' @keywords internal
create.cvg.text <- function(stats, selected.group = NULL, ident = NULL) {
    if (length(selected.group) == 0 || "all" %in% selected.group) {
        selected.group <- "Total"
    } 
    if (selected.group[1] == "Total") {
        out.selected.group <- "" # for output: change id
    } else {
        out.selected.group <- selected.group
    }
    if (length(ident) == 0) { # identifier for coverage type
        ident <- "Coverage"
    }
    idx <- which(stats$Group %in% selected.group)
    if (length(idx) != 0) {
        stats <- stats[idx, ]
        # differentiate fw & rev cvg
        text <- paste0(ident, ": ", stats$N_primer, " primers (", stats$N_primer_fw, " fw/", stats$N_primer_rev, 
            " rev) cover ", stats$Coverage, " ", out.selected.group, " template sequences.")
    } else {
        text <- "Templates not covered."
    }
    return(text)
}
#' Check for Report Dependencies.
#'
#' Checks whether the dependencies for rmarkdown::render() are fulfilled.
#'
#' @return A logical vector giving the dependency availability status.
#' @keywords internal
check_report_deps <- function() {
	pandoc.available <- Sys.which("pandoc") != ""
	latex.available <- Sys.which("pdflatex") != ""
	out <- c("Pandoc" = pandoc.available, "LateX" = latex.available)
	return(out)
}
#' Sanitiziation of Filename.
#'
#' Ensures that a filename is valid for the file system.
#'
#' @param path The path to the file to be sanitized, without file extension.
#' @param suffix The suffix (e.g. ".png") of a file.
#' @param sub_char The character used to replacing disallowed chars.
#' @return The sanitized filename
#' @keywords internal
sanitize_path <- function(path, suffix = '', sub.char = "_") {
  if (grepl('[^~:_./\\[:alnum:]-]', path)) {
    warning('Replaced special characters in path "', path, '" -> "',
            path <- gsub('[^~:_./\\[:alnum:]-]', sub.char, path), '"')
  }
  # replace . with sub.char except ../ and ./
  s = strsplit(path, '[/\\\\]')[[1L]]
  i = (s != '.') & (s != '..') & grepl('\\.', s)
  if (any(i)) {
    s[i] = gsub('\\.', sub.char, s[i])
    path = paste(s, collapse = '/')
    warning('Dots in path replaced with", sub.char, "("', path, '")')
  }
  out <- paste0(path, suffix)
  return(out)
}
#' Creation of a Filename for Reports.
#'
#' Creates the filename for reports.
#'
#' @param report.name The identifier for the report type.
#' @param sample.name The identifier of the sample that was analyzed.
#' @return A character vector.
#' @keywords internal
get_report_fname <- function(report.name, sample.name) {
    date <- format(Sys.time(), "%Y-%m-%d")
    out <- paste(report.name, "_", sample.name, "_", date, sep = "")
    out <- sanitize_path(out, ".pdf")
    return(out)
}
#' Renders an rmarkdown file using Pandoc.
#'
#' Creates a PDF report using rmarkdown and Pandoc by passing the specified
#' \code{params} to the markdown file given by \code{report_template} and 
#' storing the PDF in \code{out.file}.
#'
#' @param params A list with parameters for the R markdown parser.
#' @param report_template A character vector giving the basename of the
#' Rmarkdown template to use for report creation.
#' @param out.file The filename of the report PDF to be created.
#' @return Creates a PDF in \code{out.file} if successful.
#' @keywords internal
render_report <- function(params, report_template, out.file) {
    deps <- check_report_deps()
	missing.deps <- names(deps)[which(!deps)]
	if (length(missing.deps) != 0) {
		msg <- paste("Cannot create PDF report. Dependencies are missing:\n",
					paste(missing.deps, collapse = ","))
		warning(msg)
		return()
	}
    report.template.dir <- system.file("extdata", "report", 
                            package = "openPrimeR")
    # copy report file to temporary directory to make sure we can write
    generation.dir <- tempdir() # dir where report is generated
    dir.copy(report.template.dir, generation.dir, overwrite = TRUE)
    tempReport <- file.path(generation.dir, report_template)
    if (rmarkdown::pandoc_available()) {
        dir.create(dirname(out.file), showWarnings = FALSE, recursive = TRUE)
        out.file <- file.path(normalizePath(dirname(out.file)), basename(out.file))
        rmarkdown::render(tempReport, output_file = out.file, 
            params = params, envir = new.env(parent = globalenv()),
            quiet = TRUE)
    } else {
        msg <- paste0("Pandoc for rmarkdown is not available on your system.", 
                "Please install it first to generate a report.")
        warning(msg)
        return()
    }
}

#' Read FASTA File.
#'
#' Reads the input FASTA file.
#'
#' @param fasta.file The path to a FASTA file.
#' @param NTs The allowed set of nucleotides.
#' @return List with vectors of chars.
#' @keywords internal
my.read.fasta <- function(fasta.file, NTs) {
    # NTs: supported NTs
    if (!file.exists(fasta.file)) {
        stop(paste("The input fasta.file did not exist at: '",
            fasta.file, "'", sep = ""))
    }
    seqs <- suppressWarnings(seqinr::read.fasta(fasta.file, forceDNAtolower = TRUE))
    # replace empty chars
    for (i in seq_along(seqs)) {
        if (" " %in% seqs[[i]]) {
            new.seq <- seqs[[i]][seqs[[i]] != " "]
            seqs[[i]] <- seqinr::as.SeqFastadna(new.seq, attr(seqs[[i]], "name"), attr(seqs[[i]], 
                "Annot"))
        }
    }
    # check whether the seqs are in the NT alphabet
    sanity.check <- unlist(lapply(seqs, function(x) all(x %in% NTs)))
    idx <- which(!sanity.check)
    # message(idx)
    if (length(idx) != 0) {
        msg <- "Some sequences contained non-supported characters. Supported are the following characters: "
        sup <- paste(NTs, collapse = ", ")
        loc <- paste(sapply(idx, function(x) attr(seqs[[x]], "name")), collapse = ", ")
        pre <- "The sequences with the following headers contained non-supported characters: "
        msg <- paste(msg, sup, ". ", pre, loc, sep = "")
        my.error("FastaAlphabetError", msg)
    }
    return(seqs)
}
#' Wrapper for the ggplot2::ggsave function.
#'
#' Saves a plot using ggplot2's ggsave function.
#'
#' @param filename The filename to store the plot.
#' @param plot The ggplot object.
#' @param ... Further arguments to the ggplot2 ggsave function.
#' @return Stores \code{p} in \code{fname}.
#' @keywords internal
my_ggsave <- function(filename, plot = ggplot2::last_plot(), ...) {
   if (length(plot) != 0) {
        check <- try(suppressMessages(ggsave(filename = filename, plot = plot, ...)))
        if (class(check) == "try-error") {
            # downgrade error to warning
            warning(attr(check, "condition") )
        }
   }
   return()
}
add_cvg_to_workbook <- function(cvg.matrix, wb, start.row, start.col, is.first.entry = FALSE) {
    # wb workbook
    # start.col The column to start writing into (counting from 1)
    # position for row/column style separators:
    # c <- seq(start.col, ncol(cvg.matrix) - 1)
    # r <- seq(start.row, start.row + nrow(cvg.matrix) - 1)
    if (is.first.entry) {
        #r <- r + 1
        openxlsx::writeData(wb, 1, cvg.matrix, 
              startCol = start.col, startRow = start.row,
              rowNames = FALSE, colNames = TRUE, borders = "all",
              borderColour = "#000000")

    } else {
        openxlsx::writeData(wb, 1, cvg.matrix, 
              startCol = start.col, startRow = start.row,
              rowNames = FALSE, colNames = FALSE, borders = "all",
              borderColour = "#000000")

    }
    # styles:
    #bodyStyle <- createStyle(border="Bottom", 
                            #borderColour = "#000000",
                            #borderStyle = "thick")
    #r <- tail(r,1)
    # doesnt work:
    #addStyle(wb, sheet = 1, bodyStyle, rows = r, cols = c, gridExpand = TRUE)
    #setColWidths(wb, 1, cols=1, widths = 21) ## set column width for row names column
    return(wb)
}
#' Creation of a Coverage XLS Spreadsheet.
#'
#' Creation of an XLS spreadsheet providing an overview of the covered
#' template sequences for each primer. Each cell in the spreadsheet
#' indicates a coverage event between a primer and template using
#' color codes. Identified coverage events are indicated by green, while
#' primer-template pairs without coverage are indicated by red.
#' In case that a primer binding condition (see \code{\link{CoverageConstraints}})
#' was active when computing the coverage, the numeric value of the
#' coverage condition is annotated for each cell.
#'
#' @param primer.df A \code{Primers} object containg primers
#' with evaluated coverage.
#' @param template.df A \code{Templates} object containing the templates
#' corresponding to \code{primer.df}.
#' @param filename A character vector providing the filename for storing the XLS file.
#' @param settings A \code{DesignSettings} object providing
#' the coverage conditions that are to be shown in the spreadsheet.
#' @return Creates a spreadsheet visualizing the primer coverage
#' and stores it under \code{filename}.
#' @export
#' @examples
#' data(Ippolito)
#' filename <- tempfile("cvg_overview", fileext = ".xls")
#' # store coverage xls file to disk
#' create_coverage_xls(primer.df, template.df, filename, settings)
create_coverage_xls <- function(primer.df, template.df, filename, settings) {
    if (length(primer.df) == 0 || nrow(primer.df) == 0) {
        return(NULL)
    }
    if (!is(settings, "DesignSettings")) {
        stop("A DesignSettings object is required for settings.")
    }
    if (!is(primer.df, "Primers")) {
        stop("Please input a 'Primers' object.")
    }
    if (!is(template.df, "Templates")) {
        stop("Please input a 'Templates' object.")
    }
    if (!"primer_coverage" %in% colnames(primer.df)) {
        warning("Please compute primer coverage first.")
        return(NULL) 
    }
    cvg.constraints <- cvg_constraints(settings)
    constraints <- names(cvg.constraints)[names(cvg.constraints) %in% c("primer_efficiency", "annealing_DeltaG")]
    constraints <- constraints[constraints %in% colnames(primer.df)]
    if (length(constraints) > 1) {
        warning("Only one constraint possible at a time. Selecting the first one.")
        constraints <- constraints[1]
    }
    # initialize workbook:
    wb <- openxlsx::createWorkbook()
    openxlsx::addWorksheet(wb, primer.df$Run[1])
    # create coverage data for individual mismatch settings:
    mm <- max(unique(unlist(sapply(strsplit(c(primer.df$Nbr_of_mismatches_fw, primer.df$Nbr_of_mismatches_rev), split = ","), function(x) as.numeric(x)))))
    mm <- c(0, seq_len(mm))
    # add mismatch column to worksheet:
    mm.column <- data.frame("Allowed_Mismatches" = unlist(lapply(mm, function(x) rep(x, nrow(primer.df)))))
    openxlsx::writeData(wb, 1, mm.column , startCol = 1, 
              rowNames = FALSE)
    #openxlsx::openXL(wb) # error already here: SET_VECTOR_ELT() can only be applied to a list, not a symbol

    start.row <- 1
    for (i in seq_along(mm)) {
        # create coverage data for one mismatch setting
        allowed.mm <- mm[i]
        fw.mm <- lapply(strsplit(primer.df$Nbr_of_mismatches_fw, split = ","), as.numeric)
        allowed.idx.fw <- lapply(fw.mm, function(x) which(x <= allowed.mm))
        rev.mm <- lapply(strsplit(primer.df$Nbr_of_mismatches_fw, split = ","), as.numeric)
        allowed.idx.rev <- lapply(fw.mm, function(x) which(x <= allowed.mm))
        # select worst case inclusion
        allowed.idx <- lapply(seq_len(nrow(primer.df)), function(x) {
                            if (primer.df$Direction[x] == "fw") {
                                allowed.idx.fw[[x]]
                            } else if (primer.df$Direction[x] == "rev") {
                                allowed.idx.rev[[x]]
                            } else { # both
                                intersect(allowed.idx.fw[[x]],allowed.idx.rev[[x]])
                            }
        })
        # update covered seqs for allowed nbr of mismatches:
        cur.df <- primer.df
        cvd.idx <- covered.seqs.to.idx(primer.df$Covered_Seqs, template.df)
        cvd <- sapply(seq_along(allowed.idx), function(x) paste(template.df$Identifier[cvd.idx[[x]][allowed.idx[[x]]]], collapse = ","))
        cur.t.df <- update_template_cvg(template.df, cur.df)
        cur.df$Covered_Seqs <- cvd
        # get coverage matrix for specific number of mismatches
        if (length(constraints) == 0) {
            cvg.matrix <- t(get.coverage.matrix(cur.df, cur.t.df))
        } else {
            # change entries to cvg constraint values: efficiency / deltaG of annealing
            # re-write constraints: select only the current constraint events
            split.vals <- strsplit(primer.df[, constraints], split = ",")
            vals <- unlist(lapply(seq_along(split.vals), function(x) paste(split.vals[[x]][allowed.idx[[x]]], collapse = ",")))
            cur.df[,constraints] <- vals
            cvg.matrix <- t(get.coverage.matrix(cur.df, cur.t.df, constraints = constraints))
            # overwrite previous 0/1 coverag entries with cvg constraint values
        }
        cnames <- c("Primer", colnames(cvg.matrix))
        cvg.matrix <- cbind(Primer = rownames(cvg.matrix), 
                                data.frame(cvg.matrix, stringsAsFactors = FALSE))
        cvg.matrix$Primer <- as.character(cvg.matrix$Primer)
        colnames(cvg.matrix) <- cnames
        # add entries to notebook
        wb <- add_cvg_to_workbook(cvg.matrix, wb, start.row = start.row, start.col = 2, is.first.entry = (i == 1))
        if (i == 1) {
            # first entry in workbook
            start.row <- start.row + nrow(primer.df) + 1 # plus the header row
        } else {
            start.row <- start.row + nrow(primer.df)
        }
    }
    # conditional formatting is applied when xls is created and based on the final values ...
    # color according to cvg.matrix entries -> 1 covered, 0 not
    coveredStyle <- openxlsx::createStyle(fontColour = "#000000", bgFill = "#00b050")
    uncoveredStyle <- openxlsx::createStyle(fontColour = "#000000", bgFill = "#ff0000")
    # row and column indices for style rules:
    r <- seq_len(nrow(mm.column)) + 1
    c <- seq_len(nrow(template.df)) + 2 # add 2 for primer and mismatch cols
    if (length(constraints) == 0) {
        match.rule <- "==1"
        mismatch.rule <- "==0"
    } else {
        if (constraints == "primer_efficiency") {
            match.rule <- paste0(">=", cvg.constraints[[constraints]]["min"])
            mismatch.rule <- paste0("<", cvg.constraints[[constraints]]["min"])
        } else { # annealing
            match.rule <- paste0("<=", cvg.constraints[[constraints]]["max"])
            mismatch.rule <- paste0(">", cvg.constraints[[constraints]]["max"])
        }
    }
    openxlsx::conditionalFormatting(wb, 1, cols = c, rows = r, rule = match.rule, style = coveredStyle)
    openxlsx::conditionalFormatting(wb, 1, cols = c, rows = r, rule = mismatch.rule, style = uncoveredStyle)
    #openxlsx::openXL(wb) # open temporary version for testing!
    openxlsx::saveWorkbook(wb, filename, overwrite = TRUE) 
}

#' Creation of a Primer PDF Report.
#' 
#' Creates a PDF report for analyzed primer sets.
#'
#' @param primers To create a report for a single primer set, please provide
#' an evaluated \code{Primers} object.
#' For creating a report comparing multiple primer sets, please provide
#' a list of \code{Primers} objects.
#' @param templates If \code{primers} is a \code{Primers} object, \code{templates} should be a \code{Templates} object.
#' If \code{primers} is a list of \code{Primers} objects, \code{templates}
#' should be a list of \code{Templates} objects of the same length as \code{primers}.
#' @param out.file A character vector giving the path to the file where the report should be stored.
#' @param settings The \code{DesignSettings} object that was used for analyzing
#' the input primers.
#' @param sample.name An identifier for your analysis. By default ( 
#' \code{NULL}), the sample identifier is selected from the
#' \code{Run} column of the input templates.
#' @param used.settings A named list (with fields \code{fw} and \code{rev}) containing the relaxed settings
#' for designing forward/reverse primers. By default (\code{NULL}), 
#' the relaxed settings are not shown in the report.
#' @param ... \code{required.cvg} (optional, default: 1), the desired coverage ratio if \code{primers} is a single primer set.
#' @return Creates a PDF file summarizing the results from analyzing one or 
#' multiple sets of primers.
#' @keywords Primers
#' @export
#' @family primer functions
#' @include primers.R templates.R
#' @note
#' Creating the report requires the external programs Pandoc (http://pandoc.org)
#' and LaTeX (http://latex-project.org).
#' @examples
#' setting.xml <- system.file("extdata", "settings", 
#'                  "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
#' settings <- read_settings(setting.xml)
#' # Creation of a report for a single primer set
#' data(Ippolito)
#' out.file.single <- tempfile("evaluation_report", fileext = ".pdf")
#' create_report(primer.df, template.df, out.file.single, settings)
#' # Creation of a report for multiple primer sets.
#' data(Comparison)
#' out.file.comp <- tempfile("comparison_report", fileext = ".pdf")
#' create_report(primer.data[1:3], template.data[1:3], out.file.comp, settings)
setGeneric("create_report", 
    function(primers, templates, out.file, settings, 
             sample.name = NULL, used.settings = NULL, ...) {
        standardGeneric("create_report")
})
#' Creation of a PDF report.
#' 
#' Creates a PDF report for a set of primers.
#'
#' @param primers An evaluated \code{Primers} object.
#' @param templates A \code{Templates} object.
#' @param out.file A character vector giving the file to store the report in.
#' @param settings A \code{DesignSettings} object.
#' @param sample.name An identifier for your analysis.
#' @param used.settings A named list (with fields "fw" and "rev") containing the forward/reverse used design settings.
#' @param required.cvg The required coverage ratio.
#' @return Creates a PDF file reporting on the input primers.
#' @keywords internal
#' @note
#' Creating the report requires the external programs Pandoc (http://pandoc.org)
#' and LaTeX (http://latex-project.org).
setMethod("create_report", c("Primers", "Templates"), function(primers, templates, 
             out.file, settings, sample.name, used.settings, required.cvg = 1) {
    mode.directionality <- get.analysis.mode(primers)
    if (length(templates) == 0 || nrow(templates) == 0 ||
        !is(templates, "Templates")) {
        warning("No available templates for creating a report; returning.")
        return()
    }
    if (length(primers) == 0 || nrow(primers) == 0 || 
        !is(primers, "Primers")) {
        warning("No available primers for creating a report; returning.")
        return()
    }
    if (is.null(sample.name)) {
        sample.name <- primers$Run[1]
    }
    if (!is(settings, "DesignSettings") ) {
        stop("Settings are invalid.")
    }
    if (any(sapply(used.settings, function(x) !is(x, "DesignSettings")))) {
        stop("used.settings are invalid.")
    }
    if (!"primer_coverage" %in% colnames(primers)) {
        warning("Please compute primer coverage first.")
        return()
    }
    if (!requireNamespace("pander")) {
        warning("Cannot create the report without 'pander'.")
        return()
    }
    # create param list for comparison report
    params <- list("primers" = primers, "templates" = templates,
                "direction" = mode.directionality, "sample" = sample.name,
                "settings" = settings,
                "used_settings" = used.settings,
                "required_cvg" = required.cvg)
    render_report(params, "report.Rmd", out.file)
    return()
})
#####
#' Creation of a PDF Report for Primer Comparison.
#' 
#' Creates a PDF report for comparing multiple primers.
#'
#' @param primers A list with evaluated \code{Primers} objects.
#' @param templates A list with \code{Templates} objects.
#' @param out.file A character vector giving the file to store the report in.
#' @param settings A \code{DesignSettings} object.
#' @param sample.name An identifier for your analysis.
#' @param used.settings A named list (with fields "fw" and "rev") containing the forward/reverse used design settings.
#' @return Creates a PDF file giving a report on the comparison of the input primers.
#' @keywords internal
#' @note
#' Creating the report requires the external programs Pandoc (http://pandoc.org)
#' and LateX (http://latex-project.org).
setMethod("create_report", c("list", "list"), function(primers, templates, 
             out.file, settings, 
             sample.name = NULL, 
             used.settings = NULL) {
    mode.directionality <- sapply(primers, function(x) get.analysis.mode(x))
    # check type of primer and template data
    if (length(primers) == 0 || length(templates) == 0) {
        stop("No primer/template data was inputted.")
    }
    template.classes <- sapply(templates, function(x) class(x))
    primer.classes <- sapply(primers, function(x) class(x))
    if (any(template.classes != "Templates") || any(primer.classes != "Primers")) {
        stop("Check types of primers/templates.")
    }
    if (is.null(sample.name)) {
        # assume templates relate to the same sample
        sample.name <- templates[[1]]$Run[1]
    }
    if (!is(settings, "DesignSettings") ) {
        stop("Settings are invalid.")
    }
    # cvg available/evaluated set?
    cvg.available <- sapply(primers, function(x) 
                        "primer_coverage" %in% colnames(x))

    if (!all(cvg.available)) {
        warning("Please compute primer coverage first.")
        return()
    }
    params <- list("primers" = primers, "templates" = templates,
                "direction" = mode.directionality, "sample" = sample.name,
                "settings" = settings)
    render_report(params, "report_comparison.Rmd", out.file)
    return()
})


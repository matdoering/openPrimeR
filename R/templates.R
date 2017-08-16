##########################
# Template functionalities
#########################

#' Validates a Template Object.
#'
#' Checks whether a Templates object is valid or not.
#'
#' @param object An input data frame to be checked for being a template data frame.
#' @return \code{TRUE}, if the object is valid, \code{FALSE} otherwise.
#' @keywords internal
validate_templates <-  function(object) {
    # specify minimal set of columns that should be present in a template data frame:
    errors <- NULL
    required.fields <- list("ID" = "character", "Header" = "character", 
                  "Identifier" = c("integer", "numeric"), "Sequence_Length" = c("integer", "numeric"),
                  "Group" = "character", "Allowed_Start_fw" = c("integer", "numeric"),
                  "Allowed_End_fw" = c("integer", "numeric"), "Allowed_Start_rev" = c("integer", "numeric"),
                  "Allowed_End_rev" = c("integer", "numeric"), "Allowed_fw" = "character", "Allowed_rev" = "character",
                  "Allowed_Start_fw_initial" = c("integer", "numeric"), "Allowed_End_fw_initial" = c("integer", "numeric"), 
                  "Allowed_Start_rev_initial" = c("integer", "numeric"),
                  "Allowed_End_rev_initial" = c("integer", "numeric"), "Sequence" = "character",
                  "Run" = "character")
    if (!is(object, "data.frame")) {
        errors <- c(errors, "Input was no data frame.")
        return(FALSE)
    }
    possible.fields <- NULL # don't check for additional fields here
    check.fields <- check_setting(possible.fields, object, required.fields)
    if (is.character(check.fields)) {
        errors <- c(errors, check.fields)
    }
    if (length(errors) != 0) {
        return(errors) 
    } else {
        # ensure that 'Run' is unique
        check.run <- length(unique(object$Run)) <= 1
        if (check.run) {
            return(TRUE)
        } else {
            msg <- "The 'Run' column may only contain a single unique value."
            return(msg)
        }
    }
}

#' The Templates Class.
#'
#' The \code{Templates} class encapsulates a data frame
#' containing the sequencs of the templates, their binding regions,
#' as well as additional information (e.g. template coverage).
#'
#' In the following you can find a description
#' of the most important columns that can be found
#' in an object of class \code{Templates}. Note that
#' angle brackets in the column names 
#' indicate the existence of multiple possibilities.
#' \describe{
#' \item{\code{ID}}{The identifiers of the templates.}
#' \item{\code{Identifier}}{The internal identifiers of the templates.}
#' \item{\code{Group}}{The identifiers of the groups that the templates belong to.}
#' \item{\code{Allowed_Start_<fw|rev>}}{The start of the interval in the templates
#' where binding is allowed for forward and reverse primers, respectively.}
#' \item{\code{Allowed_End_<fw|rev>}}{The end of the interval in the templates
#' where binding is allowed for forward and reverse primers, respectively.}
#' \item{\code{Allowed_<fw|rev>}}{The template sequence where binding
#' is allowed for forward and reverse primers, respectively.}
#' \item{\code{Run}}{An identifier for the set of template sequences.}
#' \item{\code{Covered_By_Primers}}{The identifiers of primers covering the templates,
#' when the template coverage has been annotated.}
#' \item{\code{primer_coverage}}{The number of primers covering the templates,
#' when the template coverage has been annotated.}
#' }
#' @name Templates-class
#' @return A \code{Templates} object, an instance of a data frame.
#' @rdname Templates-class
#' @exportClass Templates
#' @export
#' @keywords Classes
#' @family template functions
#' @seealso \code{\link{read_templates}} for loading template sequences,
#' \code{\link{assign_binding_regions}} for adjusting the primer binding regions,
#' \code{\link{update_template_cvg}} for setting the template coverage,
#' \code{\link{plot_template_cvg}} for plotting the template coverage,
#'
#' @examples
#' fasta.file <- system.file("extdata", "IMGT_data", "templates", 
#'      "Homo_sapiens_IGH_functional_exon.fasta", package = "openPrimeR")
#' hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")
#' template.df <- read_templates(fasta.file, hdr.structure, "|", "GROUP")
setClass("Templates",
   contains = c("data.frame"), validity = validate_templates)


setMethod("initialize", "Templates",
    function(.Object, ...) {
        # just call data frame constructor
        .Object <- callNextMethod(.Object, ...)
    }
) 
setMethod("show", "Templates", function(object) {
    # overwrite the 'print' function using data.table's show
    my_show_df(asS3(object))
})
#' @name Templates
#' @rdname Templates-class
#' @export
#' @param ... A data frame fulfilling the structural requirements
#' for initializing a \code{Templates} object.
Templates <- function(...) new("Templates", ...) 

#' cbind for Template class.
#'
#' Ensures that the cbind result has the appropriate class.
#'
#' @param ... Parameters for cbind function.
#' @return Column binded Templates data frame.
#' @keywords internal
#' @export
#' @examples
#' data(Ippolito)
#' template.df <- cbind(template.df, seq_len(nrow(template.df)))
cbind.Templates <- function(...) {
    df <- cbind.data.frame(...)
    df <- Templates(df)
    return(df)
}
#' rbind for Template class.
#'
#' Ensures that the rbind result has the appropriate class.
#'
#' @param ... Parameters for rbind function.
#' @return Row-binded Templates data frame.
#' @keywords internal
#' @export
#' @examples
#' data(Ippolito)
#' template.df <- rbind(template.df, template.df)
rbind.Templates <- function(...) {
    df <- rbind.data.frame(...)
    df <- Templates(df)
    return(df)
}

#' S4 cbind for Templates.
#'
#' S4 cbind function for Templates data frame.
#'
#' @export
#' @return Cbinded template data frame.
#' @rdname Templates-method
#' @param x The Templates data frame.
#' @param y Another data frame.
#' @keywords internal
#' @examples
#' data(Ippolito)
#' template.df <- cbind2(template.df, seq_len(nrow(template.df)))
setMethod("cbind2", "Templates", 
    # need this in case of S3 dispatch fails.
    function(x, y, ...) {
        df <- cbind.data.frame(x, y, ...)
        df <- Templates(df)
        return(df)
    }
)
#' S4 rbind for Templates.
#'
#' S4 rbind function for templates.
#'
#' @export
#' @rdname Templates-method
#' @return Rbinded template data frame.
#' @keywords internal
#' @examples
#' data(Ippolito)
#' template.df <- rbind2(template.df, template.df)
setMethod("rbind2", "Templates", 
    # need this in case of S3 dispatch fails.
    function(x, y, ...) {
        df <- rbind.data.frame(x, y, ...)
        df <- Templates(df)
        return(df)
    }
)

#' Slicing Operator for Templates.
#'
#' Slicing of Templates data frame object.
#'
#' @param x The Template data frame.
#' @param i The row index.
#' @param j The column index.
#' @param ... Other arguments to the slice operator.
#' @param drop Simplify data frame?
#' @exportMethod [
#' @rdname Templates-method
#' @keywords internal
#' @aliases [,Templates-method
#' @return Subsetted template data frame.
#' @examples
#' data(Ippolito)
#' template.df <- template.df[1:2,]
setMethod("[", "Templates",
    function(x, i, j, ..., drop = TRUE) {
        if (missing(drop)) {
            df <- asS3(x)[i, j, ...]
        } else {
            df <- asS3(x)[i, j, ..., drop = drop]
        }
        if (class(df) == "data.frame") {
            # only set my class if we haven't simplified.
            #options("show.error.messages" = FALSE)
            t.df <- suppressWarnings(try(Templates(df), silent = TRUE))
            #options("show.error.messages" = TRUE)
            if (class(t.df) == "try-error") {
                # we removed some crucial columns -> No Templates object anymore.
                t.df <- asS3(df)
            }
            df <- t.df
        }
        return(df)
    }
)
#' Dollar Operator for Templates.
#'
#' Set a column in a Templates data frame.
#'
#' @exportMethod $<-
#' @rdname Templates-method
#' @param name The name of the column.
#' @param value The values of the column.
#' @return Templates with replaced column.
#' @keywords internal
#' @examples
#' data(Ippolito)
#' template.df$ID[1] <- "newID"
setMethod("$<-", "Templates", 
    function(x, name, value) {
        df <- asS3(x)
        eval(parse(text = paste("df$", name, " <- value", sep = "")))
        df <- Templates(df)
        return(df)
    }
)
#' Removal of Partial Sequences.
#'
#' Removes partial template sequences.
#'
#' @param template.df Template data frame.
#' @param Header keywords indiating templates that should be excluded.
#' @return Template data frame with partial sequences removed.
#' @keywords internal
remove.seqs.by.keyword <- function(template.df, keyword = "partial") {
    if (length(keyword) == 0) {
        # nothing to remove
        return(template.df)
    }
    rm.idx <- unique(unlist(lapply(keyword, function(x) grep(x, template.df$Header))))
    if (length(rm.idx) != 0) {
        template.df <- template.df[-rm.idx, ]
    }
    return(template.df)
}
#' Parse FASTA headers
#'
#' Parses the headers of a FASTA file.
#'
#' @param hdr The headers (> entries) of a FASTA File.
#' @param delim The delimiter used to separate distinct fields in the headers.
#' For example, \code{|} for headers such as > E.coli | GeneX | ...
#' @param hdr.str Names of the fields appearing in the header.
#' @param id.column Field in the header to be used used as an identifier for the sequences.
#' @return Data frame with structured information from the headers.
#' @keywords internal
parse.header <- function(hdr, delim, hdr.str, id.column) {
    result <- data.frame(Header = hdr, stringsAsFactors = FALSE)
    data.df <- NULL
    if (length(hdr.str) != 0 && delim != "") {
        s <- strsplit(gsub(">", "", hdr), split = delim, fixed = TRUE)
        data.df <- data.frame(matrix(rep(NA, length(hdr.str) * length(hdr)), ncol = length(hdr.str)), 
            stringsAsFactors = FALSE)
        colnames(data.df) <- hdr.str
        for (i in seq_along(s)) {
            h <- s[[i]]
            if (length(hdr.str) > length(h)) {
                my.error("TemplateHeaderStructure", "Header is too short for given header structure!")
            }
            elements <- h[seq_along(hdr.str)]
            data.df[i, ] <- elements
        }
    } else if (length(hdr.str) != 0 && delim == "") {
        # no delimiter or header structure given -> use the complete header info for the first specified field
        s <- gsub(">", "", hdr)
        data.df <- data.frame(Info = s, stringsAsFactors = FALSE)
        colnames(data.df) <- hdr.str[1]
    } 
    if (length(data.df) != 0) {
        result <- cbind(result, data.df)
    }
    if (!id.column %in% colnames(result)) {
        my.warning("TemplateIDColNotFound", paste("ID column: ", id.column, " was not found in the header of the template. Using first available header variable as ID.", 
            sep = ""))
        result <- cbind(ID = result[, 1], result, stringsAsFactors = FALSE)
    } else {
        result <- cbind(ID = result[, id.column], result, stringsAsFactors = FALSE)
    }
    # create an identifier to match to leader sequences: should be unique and be
    # present in leader file
    if (!"Group" %in% colnames(result)) {
        result$Group <- "default"
    }
    return(result)
}


#' Input of Template Sequences.
#' 
#' Read one or multiple files with template sequences in FASTA or CSV format.
#'
#' When supplying a FASTA file with template sequences,
#' the input arguments \code{hdr.structure}, \code{delim}, \code{id.column}, 
#' \code{rm.keywords}, \code{remove.duplicates}, \code{fw.region}, 
#' \code{rev.region}, \code{gap.character}, and \code{run} are utilized.
#' Most importantly, \code{hdr.structure} and \code{delim} should
#' match the FASTA header structure.
#' To learn more about setting the primer binding regions, consider the
#' \code{\link{assign_binding_regions}} function.
#' In contrast, when supplying a CSV file with template sequences,
#' the data are loaded without performing any modifications because the CSV file
#' should represent an object of class \code{\link{Templates}}, which 
#' can be stored using the \code{\link{write_templates}} function.
#'
#' @param template.file Path to one or multiple FASTA or CSV files 
#' containing the template sequences.
#' @param hdr.structure A character vector describing the information contained in the FASTA headers. In case that the headers of \code{fasta.file} contain
#' template group information, please include the keyword "GROUP" in
#' \code{hdr.structure}. If the numer of elements provided via \code{hdr.structure}
#' is shorter than the actual header structure, the missing fields are ignored.
#' @param delim Delimiter for the information in the FASTA headers.
#' @param id.column Field in the header to be used as the identifier of individual
#' template sequences.
#' @param rm.keywords A vector of keywords that are used to remove templates whose headers contain any of the keywords.
#' @param remove.duplicates Whether duplicate sequence shall be removed.
#' @param fw.region The positional interval from the template 5' end specifying the  
#' binding sites for forward primers. The default \code{fw.region} is set to
#' the first 30 bases of the templates.
#' @param rev.region The positional interval from the template 3' end specifying
#' the binding sites for reverse primers. The default \code{rev.region}
#' is set to the last 30 bases of the templates.
#' @param gap.character The character in the input file representing gaps. 
#' Gaps are automatically removed upon input and the default character is "-".
#' @param run An identifier for the set of template sequences. By default,
#' \code{run} is \code{NULL} and its value is set via \code{template.file}.
#' @return An object of class \code{Templates}.
#' @export
#' @family template functions
#' @keywords Templates
#' @examples
#' # Load templates from a FASTA file
#' fasta.file <- system.file("extdata", "IMGT_data", "templates", 
#'           "Homo_sapiens_IGH_functional_exon.fasta", package = "openPrimeR")
#' hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")
#' template.df.fasta <- read_templates(fasta.file, hdr.structure, "|", "GROUP")
#' # Load mutliple FASTA files
#' fasta.files <- c(fasta.file, fasta.file)
#' template.df.fastas <- read_templates(fasta.files, hdr.structure, "|", "GROUP")
#' # Load templates from a previously stored CSV file
#' csv.file <- system.file("extdata", "IMGT_data", "comparison", 
#'                "templates", "IGH_templates.csv", package = "openPrimeR")
#' template.df.csv <- read_templates(csv.file)
#' # Load multiple CSV files:
#' csv.files <- c(csv.file, csv.file)
#' template.df.csvs <- read_templates(csv.files)
#' # Load a mixture of FASTA/CSV files:
#' mixed.files <- c(csv.file, fasta.file)
#' template.data <- read_templates(mixed.files)
read_templates <- function(template.file, hdr.structure = NULL, delim = NULL, id.column = NULL, 
                           rm.keywords = NULL, remove.duplicates = FALSE, fw.region = c(1,30), 
                           rev.region = c(1,30), gap.character = "-", run = NULL) {

    if (length(template.file) > 1) {
        template.df <- read_templates_multiple(template.file, 
            hdr.structure = hdr.structure, delim = delim, 
            id.column = id.column, rm.keywords = rm.keywords,
            remove.duplicates = remove.duplicates, 
            fw.region = fw.region, rev.region = rev.region, 
            gap.character = gap.character, run = run)
    } else {
        template.df <- read_templates_single(template.file, 
            hdr.structure = hdr.structure, delim = delim, 
            id.column = id.column, rm.keywords = rm.keywords,
            remove.duplicates = remove.duplicates, 
            fw.region = fw.region, rev.region = rev.region, 
            gap.character = gap.character, run = run)

    }
    return(template.df)
}
#' Read Template CSV File
#'
#' Reads an input template CSV file.
#'
#' @param fname The filename of the input template CSV file.
#'
#' @return A template data frame.
#' @keywords internal
read_templates_csv <- function(fname) {
    if (!file.exists(fname)) {
        stop("Templates CSV file at '", fname, "' could not be found.")
    }
    templates <- try(read.csv(fname, stringsAsFactors = FALSE, row.names = NULL), silent = TRUE)
    if (class(templates) == "try-error") {
        my.error("TemplateFormatIncorrect", 
            paste0("Could not read the input file: ",
            " '", fname, "' as a CSV. Please check your input!"))
    }
    char.cols <- c("Allowed_fw", "Allowed_rev", 
                    "Covered_By_Primers", "Covered_By_Primers_fw",
                    "Covered_By_Primers_rev")
    num.cols <- c("Allowed_Start_fw","Allowed_Start_rev", 
                    "Allowed_End_fw", "Allowed_End_rev",
                    "Allowed_Start_fw_initial", "Allowed_Start_rev_initial",
                    "Allowed_End_fw_initial", "Allowed_End_rev_initial",
                    "primer_coverage", "primer_coverage_fw", 
                    "primer_coverage_rev")
    c.cols <- char.cols[char.cols %in% colnames(templates)]
    templates[, c.cols] <- apply(templates[, c.cols], 2, as.character)
    n.cols <- num.cols[num.cols %in% colnames(templates)]
    templates[, n.cols] <- apply(templates[,n.cols], 2, as.numeric)
    # set NA to empty string for character columns since "" is converted to NA by write.csv()
    templates[, c.cols][is.na(templates[, c.cols])] <- c("") 
    templates <- try(Templates(templates))
    if (class(templates) == "try-error") {
        my.error("TemplateFormatIncorrect", 
            paste0("Could not construct a 'Templates' object from the",
                " CSV input file '", fname, "'. Please check that",
                " the structure of the input file is correct."), silent = TRUE)
    }
    return(templates)
}

#' Input of a Single Template File.
#' 
#' Read template sequences from a FASTA or CSV file.
#'
#' When supplying a FASTA file with template sequences,
#' the input arguments \code{hdr.structure}, \code{delim}, \code{id.column}, 
#' \code{rm.keywords}, \code{remove.duplicates}, \code{fw.region}, 
#' \code{rev.region}, \code{gap.character}, and \code{run} are utilized.
#' Most importantly, \code{hdr.structure} and \code{delim} should
#' match the FASTA header structure.
#' When supplying a CSV file with template sequences,
#' the data are loaded without any modification and the CSV file
#' should represent an object of class \code{\link{Templates}}, which 
#' can be stored using the \code{\link{write_templates}} function.
#'
#' @param template.file Path to a FASTA or CSV file containing the template sequences.
#' @param hdr.structure A character vector describing the information contained in the FASTA headers. In case that the headers of \code{fasta.file} contain
#' template group information, please include the keyword "GROUP" in
#' \code{hdr.structure}.
#' @param delim Delimiter for the information in the FASTA headers.
#' @param id.column Field in the header to be used as the identifier.
#' @param rm.keywords A vector of keywords that are used to remove templates whose headers contain any of the keywords.
#' @param remove.duplicates Whether duplicate sequence shall be removed.
#' @param fw.region The positional interval from the template 5' end specifying the 
#' binding sites for forward primers.
#' @param rev.region The positional interval from the template 3' end specifying
#' the binding sites for reverse primers.
#' @param gap.character The character in the input file representing gaps.
#' Gaps are automatically removed upon input.
#' @param run An identifier for the template sequences.
#' @return An object of class \code{Templates}.
#' @keywords internal
read_templates_single <- function(template.file, hdr.structure = NULL, delim = NULL, id.column = NULL, 
                           rm.keywords = NULL, remove.duplicates = FALSE, fw.region = c(1,30), 
                           rev.region = c(1,30), gap.character = "-", run = NULL) {
    template.df <- try(read_templates_fasta(template.file, hdr.structure, 
                       delim, id.column, rm.keywords, remove.duplicates, 
                       fw.region, rev.region, gap.character, run), silent = TRUE)
    if (class(template.df) == "try-error") {
        template.df <- try(read_templates_csv(template.file), silent = TRUE)
        if (class(template.df) == "try-error") {
            my.error("TemplateFormatIncorrect", 
                paste0("Unsupported template input file type or error reading data for file: '", template.file, "'"))
        }    
    }
    # check template quality: degeneration
    degen <- score_degen(strsplit(template.df$Sequence, split = ""), gap.char = gap.character)
    degen.idx <- which(degen > 2^10)
    if (length(degen.idx) > 0) {
        msg <- paste("At least one template sequence was highly degenerate.",
                    "Please check the quality of the input template sequences.",
                    "The following sequences will not be disambiguated:",
                    paste(template.df$ID[degen.idx], collapse = ","))
        my.warning("TemplateDegeneration", msg)
    }
    # ensure that IDs of templates are unique
    dup.idx <- which(duplicated(template.df$ID))
    if (length(dup.idx) != 0) {
        # make IDs unique
        non.unique.IDs <- unique(template.df$ID[dup.idx])
        for (ID in non.unique.IDs) {
            idx <- which(template.df$ID == ID)
            new.IDs <- paste0(ID,  "_", seq_along(idx))
            template.df$ID[idx] <- new.IDs
        }
    }
    return(template.df)
}
####
#' Input of Multiple Templates.
#'
#' Reads multiple template CSV/FASTA files.
#'
#' @param filenames Names of FASTA/CSV files containing template data.
#' @param hdr.structure A character vector describing the information contained in the FASTA headers. In case that the headers of \code{fasta.file} contain
#' template group information, please include the keyword "GROUP" in
#' \code{hdr.structure}.
#' @param delim Delimiter for the information in the FASTA headers.
#' @param id.column Field in the header to be used as the identifier.
#' @param rm.keywords A vector of keywords that are used to remove templates whose headers contain any of the keywords.
#' @param remove.duplicates Whether duplicate sequence shall be removed.
#' @param fw.region The positional interval from the template 5' end specifying the 
#' binding sites for forward primers.
#' @param rev.region The positional interval from the template 3' end specifying
#' the binding sites for reverse primers.
#' @param gap.character The character in the input file representing gaps.
#' Gaps are automatically removed upon input.
#' @param run An identifier for the template sequences.
#' @return A list containing objects of class \code{Templates}.
#' @keywords internal
read_templates_multiple <- function(filenames,  hdr.structure = NULL, delim = NULL, id.column = NULL, 
                           rm.keywords = NULL, remove.duplicates = FALSE, fw.region = c(1,30), 
                           rev.region = c(1,30), gap.character = "-", run = NULL) {
    data <- lapply(filenames, function(x) read_templates_single(x,
            hdr.structure = hdr.structure, delim = delim, 
            id.column = id.column, rm.keywords = rm.keywords,
            remove.duplicates = remove.duplicates, 
            fw.region = fw.region, rev.region = rev.region, 
            gap.character = gap.character, run = run))
    data <- set.run.names(data, filenames)
    return(data)
}

#' Input of Template Sequences.
#' 
#' Read template sequences from a FASTA file.
#'
#' @param fasta.file Path to a FASTA file containing the template sequences.
#' @param hdr.structure Names describing the information contained in the FASTA headers. In case that the headers of \code{fasta.file} contain
#' template group information, please include the keyword "GROUP" in
#' \code{hdr.structure}.
#' @param delim Delimiter for the information in the FASTA headers.
#' @param id.column Field in the header to be used as the identifier.
#' @param rm.keywords A vector of keywords that are used to remove templates whose headers contain any of the keywords.
#' @param remove.duplicates Whether duplicate sequence shall be removed.
#' @param fw.region The positional interval from the template 5' end specifying the 
#' binding sites for forward primers.
#' @param rev.region The positional interval from the template 3' end specifying
#' the binding sites for reverse primers.
#' @param gap.character The character in the input file representing gaps.
#' Gaps are automatically removed upon input.
#' @param run An identifier for the template sequences.
#' @return An object of class \code{Templates}.
#' @keywords internal
#' @examples
#' fasta.file <- system.file("extdata", "IMGT_data", "templates", 
#'               "Homo_sapiens_IGH_functional_exon.fasta", package = "openPrimeR")
#' hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")
#' template.df <- read_templates(fasta.file, hdr.structure, "|", "GROUP")
read_templates_fasta <- function(fasta.file, hdr.structure = NULL, delim = NULL, id.column = NULL, 
                           rm.keywords = NULL, remove.duplicates = FALSE, fw.region = c(1,30), 
                           rev.region = c(1,30), gap.character = "-", run = NULL) {
    allowed.nts <- c(tolower(names(Biostrings::IUPAC_CODE_MAP)), gap.character)
    seqs <- my.read.fasta(fasta.file, allowed.nts)
    if (is.null(seqs)) {
        return(NULL)
    }
    if (is.null(run)) {
        # use filename w/o extension as identifier
        run <- sub("^([^.]*).*", "\\1", basename(fasta.file))
    }
    if (length(id.column) != 0) { 
        if (length(id.column) != 1) {
            stop("Please specify only a single id.column.")
        }
        id.col.ok <- any(grepl(id.column[1], hdr.structure))
        if (!id.col.ok) {
            my.error("ID_Column_Not_Found", paste0("The specified",
                "ID column  was not part of the header structure."))
        }
    } else {
        # use the first hdr.structure field as the ID if not id.col is specified
        id.column <- hdr.structure[1]
    }
    if (length(hdr.structure) != 0 && length(delim) == 0) {
        warning("Header structure provided, but no delimiter. Cannot read meta information!")
    }
    if (length(hdr.structure) == 0 && length(delim) != 0) {
        warning("Header delimiter supplied, but no header structure. Cannot read meta information!")
    }
    if (length(delim) == 0) {
        delim <- ""
    }
    if (length(hdr.structure) == 0) {
        # no hdr.structure supplied
        id.column <- "HEADER"  # use the full header as the ID column
    }
    # retrieve seqs/replace gap char
    s <- sapply(seqs, function(x) paste(x, collapse = ""))
    hdr <- sapply(seqs, function(x) attr(x, "Annot"))
    # determine IDs and groups parse hdr structure
    hdr.str <- paste(substr(hdr.structure, 1, 1), 
                substr(tolower(hdr.structure), 2, nchar(hdr.structure)), 
                        sep = "")
    id.column <- paste(substr(id.column, 1, 1), substr(tolower(id.column), 2, nchar(id.column)), 
        sep = "")  # rename according to renaming in parse.hdr.structure
    ids <- parse.header(hdr, delim, hdr.str, id.column)  # data frame from parsed header
    if (length(ids) == 0) {
        return(NULL)
    }
    d <- data.frame(ids, Identifier = 1:length(s), InputSequence = s, Sequence_Length = nchar(s), 
        stringsAsFactors = FALSE)
    rownames(d) <- NULL
    # annotate gene groups: try to parse IMGT header structure
    gene.groups <- try(parse.IMGT.gene.groups(d$Group))
    if (class(gene.groups[1]) == "try-error") {
            gene.groups <- NULL
    }
    if (!is.null(gene.groups)) {
        # IMGT header structure -> get some more info
        groups <- gene.groups$Gene_Group 
        if (!all(is.na(groups))) {
            d$Group <- groups
        }
    }
    # exclude partial seqs
    if (length(rm.keywords) != 0) {
        d <- remove.seqs.by.keyword(d, rm.keywords)
    }
    # exclude duplicates
    if (remove.duplicates) {
        idx <- which(duplicated(d$InputSequence))
        if (length(idx) != 0) {
            d <- d[-idx,]
        }
    }
    # add default binding regions:
    d$Sequence <- ungap_sequence(d$InputSequence, gap.character)
    default.leaders <- create.uniform.leaders(fw.region, rev.region, d, gap.character)
    d <- cbind(d, default.leaders)
    d$Run <- run
    # re-order:
    last.cols <- c("Sequence", "InputSequence", "Run")
    cols <- setdiff(colnames(d), last.cols)
    cols <- c(cols, last.cols)
    d <- d[, cols]
    template.df <- Templates(d)
    return(template.df)
}
#' Read Sequences.
#'
#' Reads an input FASTA file.
#'
#' @param fasta.file The path to a FASTA file.
#' @param The character indicating gaps in sequences.
#' @return  A data frame with sequences.
#' @keywords internal
read.sequences <- function(fasta.file, gap.character) {
    if (length(fasta.file) == 0) {
        return(NULL)
    }
    # supress warning when there's an extra-newline tryCatch(
    allowed.nts <- c(tolower(names(Biostrings::IUPAC_CODE_MAP)), gap.character)
    seqs <- my.read.fasta(fasta.file, allowed.nts)
    s <- sapply(seqs, function(x) paste(x, collapse = ""))
    hdr <- sapply(seqs, function(x) attr(x, "Annot"))
    d <- data.frame(Identifier = 1:length(s), Header = hdr, Sequence = s, Sequence_Length = nchar(s), 
        stringsAsFactors = FALSE)
    rownames(d) <- NULL
    return(d)
}
#' Annotation of Template Coverage.
#'
#' Annotates the templates with coverage information.
#'
#' @param template.df An object of class \code{Templates}.
#' @param primer.df An object of class \code{Primers} containing
#' primers with annotated coverage that are to be used to update 
#' the template coverage in \code{template.df}.
#' @param mode.directionality Directionality of primers.
#' The default is \code{NULL}, which means that the directionality
#' of primers is identified automatically.
#' @return An object of class \code{Templates} with updated coverage columns.
#' @export
#' @family template functions
#' @examples
#' data(Ippolito)
#' template.df <- update_template_cvg(template.df, primer.df)
update_template_cvg <- function(template.df, primer.df, 
                                mode.directionality = NULL) {

    if (length(template.df) == 0 || length(primer.df) == 0) {
        # nothing to update ...
        return(NULL)
    }
    if (!is(primer.df, "Primers")) {
        stop("Please supply a valid primer data frame.")
    }
    if (!is(template.df, "Templates")) {
        stop("Please provide a valid template data frame.")
    }
    # check whether primers and templates correspond to one another
    if (!check_correspondence(primer.df, template.df)) {
        msg <- paste("Could not match all coverage events to the input template data frame.",
            "Please check whether the primers and templates correspond.")
        stop(msg)
    }
    if (is.null(mode.directionality)) {
        mode.directionality <- get.analysis.mode(primer.df)
    } else {
        mode.directionality <- match.arg(mode.directionality, c("fw", "rev", "both"))
    }
    cvg.matrix <- get.coverage.matrix(primer.df, template.df)  # template x primers: 1/0
    if (length(cvg.matrix) != 0) {
        # require any primer to cover a template
        # for direction "both": covered seqs are fw AND rev primers 
        idx <- lapply(seq_len(nrow(cvg.matrix)), function(x) 
            which(cvg.matrix[x,] != 0))  # idx of covering primers per template
    } else {
        idx <- NULL
    }
    # add columns:
    cols <- c("Covered_By_Primers")
    extra.cols <- c(paste(cols, "_fw", sep = ""), paste(cols, "_rev", sep = ""))
    all.cols <- c(cols, extra.cols)
    for (i in seq_along(all.cols)) {
        template.df[, all.cols[i]] <- rep("", nrow(template.df))
    }
    cols <- c("primer_coverage")
    extra.cols <- c(paste(cols, "_fw", sep = ""), paste(cols, "_rev", sep = ""))
    all.cols <- c(cols, extra.cols)
    for (i in seq_along(all.cols)) {
        template.df[, all.cols[i]] <- rep(0, nrow(template.df))
    }
    if (length(idx) != 0) {
        fw.primer.idx <- which(primer.df$Forward != "")
        rev.primer.idx <- which(primer.df$Reverse != "")
        if (mode.directionality == "both") { # require fw & rev binding events
            # coverage by both a foward and a reverse primer?
            both.idx <- sapply(idx, function(x) any(x %in% fw.primer.idx) & any(x %in% rev.primer.idx))
            primer.IDs <- lapply(seq_along(idx), function(x) if(both.idx[x]) primer.df$Identifier[idx[[x]]] else numeric(0))
        } else {
            primer.IDs <- lapply(idx, function(x) primer.df$Identifier[x])
            
        }
        template.df$Covered_By_Primers <- sapply(primer.IDs, function(x) paste(x, collapse = ","))
        # n.b.: for primers with two directions, template cvg is total nbr of fw/rev primers!
        template.df$primer_coverage <- sapply(primer.IDs, length) 
        # individual direction coverage:
        new.entries <- data.frame(
                Covered_By_Primers_fw = "", 
                primer_coverage_fw = 0, 
                Covered_By_Primers_rev = "", 
                primer_coverage_rev = 0, stringsAsFactors = FALSE)
        col.idx <- which(!colnames(new.entries) %in% colnames(template.df))
        if (length(col.idx) != 0) {
            # add missing fields 
            template.df <- cbind(template.df, new.entries[,col.idx])
        }
        if (mode.directionality == "fw" || mode.directionality == "both") {
            idx.fw <- lapply(idx, function(x) intersect(x, fw.primer.idx))
            primer.IDs.fw <- lapply(idx.fw, function(x) primer.df$Identifier[x])
            template.df$Covered_By_Primers_fw <- sapply(primer.IDs.fw, function(x) paste(x, 
                collapse = ","))
            template.df$primer_coverage_fw <- sapply(primer.IDs.fw, length)
        }
        if (mode.directionality == "rev" || mode.directionality == "both") {
            idx.rev <- lapply(idx, function(x) intersect(x, rev.primer.idx))
            primer.IDs.rev <- lapply(idx.rev, function(x) primer.df$Identifier[x])
            template.df$Covered_By_Primers_rev <- sapply(primer.IDs.rev, function(x) paste(x, 
                collapse = ","))
            template.df$primer_coverage_rev <- sapply(primer.IDs.rev, length)
        }
    } 
    # don't order by cvg -> is bad for all indecing steps that depend on uniform order
    #o <- order(template.df$primer_coverage, decreasing = TRUE)
    #template.df <- template.df[o, ]
    return(template.df)
}

#' Adjustment of Existing Binding Regions.
#'
#' Adjusts the existing annotation of binding regions by specifying
#' a new binding interval relative to the existing binding region. 
#'
#' The new binding intervals provided by \code{fw} and \code{rev}
#' for forward and reverse primers, respectively, are provided
#' relative to the existing definition of binding regions in \code{template.df},
#' which can be set using \code{\link{assign_binding_regions}}.
#' When specifying relative positions, position \code{0} is defined as the first position after the end of the existing
#' binding region. Hence, negative positions relate to regions within the
#' existing binding region while non-negative values 
#' extend the binding region further.
#'
#' @param template.df A \code{Templates} object providing the template sequences
#' for which the binding regions shall be adjusted.
#' @param region.fw Interval of new binding regions relative to the forward binding region defined in \code{template.df}.
#' @param region.rev Interval of new binding regions relative to the reverse binding region defined in \code{template.df}
#' @return A \code{Templates} object with updated binding regions.
#' @export
#' @keywords Templates
#' @examples
#' data(Ippolito)
#' # Extend the binding region by one position
#' relative.interval <- c(-max(template.df$Allowed_End_fw), 0)
#' template.df.adj <- adjust_binding_regions(template.df, relative.interval)
#' # compare old and new annotations:
#' head(cbind(template.df$Allowed_Start_fw, template.df$Allowed_End_fw))
#' head(cbind(template.df.adj$Allowed_Start_fw, template.df.adj$Allowed_End_fw))
adjust_binding_regions <- function(template.df, region.fw, region.rev) {
    if (!is(template.df, "Templates")) {
        stop("'template.df' should be a 'Templates' object.")
    }
    if (missing(region.fw)) {
        region.fw <- NULL
    }
    if (missing(region.rev)) {
        region.rev <- NULL
    }
    cols <- c("Allowed_Start_fw", "Allowed_End_fw", "Allowed_Start_rev", "Allowed_End_rev")
    if (length(region.fw) != 0) {
        if (length(region.fw) != 2) {
            stop("'region.fw' wasn't an interval.")
        }
        min <- region.fw[1]
        max <- region.fw[2]
        template.df <- update.individual.binding.region(min, max, template.df, "fw")
    }
    if (length(region.rev) != 0) {
        if (length(region.rev) != 2) {
            stop("'region.rev' wasn't an interval.")
        }
        min <- region.rev[1]
        max <- region.rev[2]
        template.df <- update.individual.binding.region(min, max, template.df, "rev")
    }
    return(template.df)
}
#' Adjustment of Existing Binding Regions for one Direction.
#'
#' Adjusts the existing annotation of binding regions by specifying
#' an interval relative to the existing binding region. 
#'
#' Position \code{0} indicates the first position after the existing
#' binding region. Hence, negative positions adjust the binding region towards
#' the existing binding regions and non-negative positionis extend
#' the existing binding region definition aways from the existing target region.
#'
#' @param min Position where binding should start.
#' @param max End position of binding.
#' @param Template data frame.
#' @param mode.directionality Directionality of primers.
#' @return Template data frame with updated binding regions.
#' @keywords internal
update.individual.binding.region <- function(min, max, template.df, mode.directionality) {
    # min/max: new setting for binding position relative to target region
    allowed.start.initial <- paste("Allowed_Start_", mode.directionality, "_initial", 
        sep = "")
    allowed.start <- paste("Allowed_Start_", mode.directionality, sep = "")
    allowed.end.initial <- paste("Allowed_End_", mode.directionality, "_initial", 
        sep = "")
    allowed.end <- paste("Allowed_End_", mode.directionality, sep = "")
    allowed.region <- paste("Allowed_", mode.directionality, sep = "")
    if (mode.directionality == "fw") {
        target.pos <- template.df[, allowed.end.initial] + 1  # absolute target position
    } else {
        target.pos <- template.df[, allowed.start.initial] - 1
    }
    if (length(min) != 0 && length(max) != 0) {
        # start
        if (mode.directionality == "fw") {
            new.start <- sapply(target.pos + min, function(x) max(x, 1))
            sel <- new.start >= nchar(template.df$Sequence)
            new.start[sel] <- nchar(template.df$Sequence)[sel]
            template.df[, allowed.start] <- new.start
        } else {
            new.end <- sapply(seq_len(nrow(template.df)), function(x) min(target.pos[x] - 
                min, nchar(template.df$Sequence[x])))
            new.end[new.end <= 0] <- 1
            template.df[, allowed.end] <- new.end
        }
        # end
        if (mode.directionality == "fw") {
            new.end <- sapply(seq_len(nrow(template.df)), function(x) min(target.pos[x] + 
                max, nchar(template.df$Sequence[x])))
            new.end[new.end <= 0] <- 1
            template.df[, allowed.end] <- new.end
        } else {
            new.start <- sapply(target.pos - max, function(x) max(x, 1))
            template.df[, allowed.start] <- new.start
        }
        # update allowed region
        template.df[, allowed.region] <- substring(template.df$Sequence, template.df[, 
            allowed.start], template.df[, allowed.end])
    }
    return(template.df)
}
#' Assignment of Template Binding Regions.
#'
#' Assigns the primer target binding regions to a set of template sequences.
#'
#' The arguments \code{fw} and \code{rev} provide data describing
#' the binding regions of the forward and reverse primers, respectively.
#' To specify binding regions for each template
#' individually, \code{fw} and \code{rev} should provide the paths to FASTA files. The headers of these
#' FASTA file should match the headers of the loaded \code{template.df} 
#' and the sequences in the files specified by \code{fw} and \code{rev} should indicate
#' the target binding regions. 
#'
#' To specify uniform binding regions,
#' \code{fw} and \code{rev} should be numeric intervals indicating the allowed
#' binding range for primers in the templates. The \code{fw} interval is 
#' specified with respect to the 5' end, while the \code{rev} interval is
#' specified with respect to the template 3' end.
#' If \code{optimize.region} is \code{TRUE}, the input binding region is
#' adjusted such that regions forming secondary structures are avoided.
#'
#' @param template.df A \code{Templates} object containing the sequences
#' for which primer binding regions should be annotated.
#' @param fw Binding regions for forward primers. Either a numeric interval indicating a uniform
#' binding range relative to the template 5' end or a path to a FASTA file providing
#' binding sequences for every template. If \code{fw} is missing, only
#' \code{rev} is considered.
#' @param rev Binding regions for reverse primers. Either a numeric interval indicating a uniform
#' binding range relative to the template 3' end or the path to a FASTA file providing
#' binding sequences for every template. If \code{rev} is missing,
#' only \code{fw} is considered.
#' @param optimize.region If \code{TRUE}, the binding regions
#' specified via \code{fw} and \code{rev} are 
#' adjusted such that binding regions that may form secondary structures are 
#' avoided. This feature requires ViennaRNA (see notes). If \code{FALSE}
#' (the default), the input binding regions are not modified.
#' @param primer.length A numeric scalar providing the probe length that is used for
#' for adjusting the primer binding regions when \code{optimize.region} is \code{TRUE}.
#' @param gap.char The character indicating gaps in aligned sequences.
#' The default is "-".
#' @return An object of class \code{Templates} with assigned binding regions.
#' @export
#' @keywords Templates
#' @family template functions
#' @note The program ViennaRNA (https://www.tbi.univie.ac.at/RNA/) must be
#' installed to perform the computations that are triggered
#' when \code{optimize.region} is set to \code{TRUE}.
#' @examples
#' data(Ippolito)
#' # Assignment of individual binding regions
#' l.fasta.file <- system.file("extdata", "IMGT_data", "templates", 
#'      "Homo_sapiens_IGH_functional_leader.fasta", package = "openPrimeR")
#' template.df.individual <- assign_binding_regions(template.df, l.fasta.file, NULL)
#' # Assign the first/last 30 bases as forward/reverse binding regions
#' template.df.uniform <- assign_binding_regions(template.df, c(1,30), c(1,30))
#' # Optimization of binding regions (requires ViennaRNA)
#' \dontrun{template.df.opti <- assign_binding_regions(template.df, c(1,30), c(1,30),
#'                      optimize.region = TRUE, primer.length = 20)}
assign_binding_regions <- function(template.df, fw = NULL, rev = NULL, optimize.region = FALSE, 
                                   primer.length = 20, gap.char = "-") {
    if (is.null(fw) && !is.null(rev)) {
        fw <- eval(parse(text = paste(class(rev), "(0)", sep = "")))
    } else if (is.null(rev) && !is.null(fw)) {
        rev <- eval(parse(text = paste(class(fw), "(0)", sep = "")))
    } else if (is.null(fw) && is.null(rev)) {
        warning("'fw' and 'rev' were both NULL. No binding regions were annotated.")
        return(template.df)
    }
    if (!is(template.df, "Templates")) {
        stop("Please supply a valid template data frame.")
    }
    if (optimize.region && length(primer.length) != 1) {
        stop("Please specify a numeric, scalar primer length.")
    }
    UseMethod("assign_binding_regions", fw)
}

#' Numeric Assignment of Binding Regions.
#'
#' Numeric S3 generic case for assigning binding regions.
#'
#' @param template.df Template data frame.
#' @param fw Binding region data forward primers.
#' @param rev Binding region data for reverse primers.
#' @param optimize.region Should the primer binding region be optimized using secondary structure prediction?
#' @param primer.length Probe length for optimizing template secondary structure.
#' @return The template data frame with assigned binding regions.
#' @export
#' @keywords internal
assign_binding_regions.numeric <- function(template.df, fw, rev, optimize.region = FALSE, 
                                           primer.length = 20, gap.char = "-") {

    if (missing(rev)) {
        rev <- NULL
    }
    if (length(fw) > 2 || length(rev) > 2 ||
        length(fw) == 1 || length(rev) == 1) {
        stop("fw or rev were no intervals.")
    }
    template.df$Sequence <- ungap_sequence(template.df$InputSequence, gap.char)
    l.seq <- create.uniform.leaders(fw, rev, template.df, gap.char)
    template.df <- add.uniform.leaders.to.seqs(template.df, l.seq)
    # optimization of regions:
    if (optimize.region) {
        if (length(fw) != 0 && length(rev) != 0) {
            mode.directionality <- "both"
        } else if (length(fw) != 0) {
            mode.directionality <- "fw"
        } else {
            mode.directionality <- "rev"
        }
        opti.regions <- optimize.template.binding.regions.dir(template.df, NULL, primer.length, mode.directionality)
        opti.intervals <- opti.regions$Intervals
        template.df <- update.binding.regions(template.df, opti.intervals)
    }
    return(template.df)
}
#' Character Assignment of Binding Regions
#'
#' Generic method for assigning the binding region using individual binding regions.
#'
#' @param template.df Template data frame.
#' @param fw FASTA file specifying the forward binding regions.
#' @param rev FASTA file specifying the reverse binding regions.
#' @param optimize.region Should the primer binding region be optimized using secondary structure prediction?
#' @param primer.length Probe length for optimizing template secondary structure.
#' @param gap.char The gap character for aligned sequences.
#' @return Template data frame with assigned binding regions.
#' @export
#' @keywords internal
assign_binding_regions.character <- function(template.df, fw, rev, optimize.region = FALSE, 
                                    primer.length = 20, gap.char = "-") {

    fw.df <- read.leaders(fw, "fw", gap.character = gap.char)
    rev.df <- read.leaders(rev, "rev", gap.character = gap.char)
    uni.leaders <- unify.leaders(fw.df, rev.df, template.df, gap.char)
    template.df <- get.leader.exon.regions(template.df, uni.leaders)
    # optimization of regions:
    if (optimize.region) {
        if (length(fw) != 0 && length(rev) != 0) {
            mode.directionality <- "both"
        } else if (length(fw) != 0) {
            mode.directionality <- "fw"
        } else {
            mode.directionality <- "rev"
        }
        opti.regions <- optimize.template.binding.regions.dir(template.df, NULL, primer.length, mode.directionality)
        opti.intervals <- opti.regions$Intervals
        template.df <- update.binding.regions(template.df, opti.intervals)
    }
    return(template.df)
}
#' Update of Binding Regions.
#'
#' Updates the binding regions in the templates by providing
#' new intervals for forward and reverse binding regions.
#'
#' @param template.df An object of class \code{Templates}.
#' @param opti.regions List with new binding intervals.
#' The list can contain the components \code{fw} and \code{rev} providing
#' numeric vectors of length 2 providing the start and end
#' of the binding regions in the templates, for forward
#' and reverse binding regions, respectively.
#' @return A \code{Templates} object with updated binding regions.
#' @keywords internal
update.binding.regions <- function(template.df, opti.regions) {
    # updates binding regions
    if ("fw" %in% names(opti.regions)) {
        if (length(opti.regions[["fw"]]) != 0) {
            starts <- sapply(opti.regions[["fw"]], function(x) x[1])
            sel <- !is.na(starts)  # only modify those where we could find a possible region ...
            template.df$Allowed_Start_fw[sel] <- starts[sel]
            template.df$Allowed_End_fw[sel] <- sapply(opti.regions[["fw"]], function(x) x[2])[sel]
            template.df$Allowed_fw[sel] <- substring(template.df$Sequence[sel], template.df$Allowed_Start_fw[sel], 
                template.df$Allowed_End_fw[sel])
        }
    }
    if ("rev" %in% names(opti.regions)) {
        if (length(opti.regions[["rev"]]) != 0) {
            starts <- sapply(opti.regions[["rev"]], function(x) x[1])
            sel <- !is.na(starts)
            template.df$Allowed_Start_rev[sel] <- starts[sel]
            template.df$Allowed_End_rev[sel] <- sapply(opti.regions[["rev"]], function(x) x[2])[sel]
            template.df$Allowed_rev[sel] <- substring(template.df$Sequence[sel], template.df$Allowed_Start_rev[sel], 
                template.df$Allowed_End_rev[sel])
        }
    }
    return(template.df)
}
#' Ungapping of Sequences.
#'
#' Removes gaps from the input sequences.
#'
#' @param seqs The input character vector with sequences
#' @param gap.char The character used to represent gaps.
#' @return \code{seqs} with gaps removed.
#' @keywords internal
ungap_sequence <- function(seqs, gap.char = "-") {
    seqs <- gsub(gap.char, "", seqs)
    return(seqs)
}
binding.interval.ungapped <- function(template.df, start, end, gap.char) {
    # adjust the binding positions with gaps to those without gaps 
    # gap count from start of seq until the end of the allowed region
    gap.count.before <- sapply(strsplit(substring(template.df$InputSequence, 1, start), split = ""),
                    function(x) length(which(x == gap.char)))
    gap.count.within <- sapply(strsplit(substring(template.df$InputSequence, start, end), split = ""),
                    function(x) length(which(x == gap.char)))
    start <- start - gap.count.before
    end <- end - gap.count.before - gap.count.within
    return(list(start, end))
}
#' Retrieval of Binding Regions
#'
#' Retrives information about individual binding regions.
#'
#' @param template.df Template data frame.
#' @param direction Identify forward and reverse.
#' @param start Start positions.
#' @param end End positions.
#' @param gap.char The character for gaps in alignments.
#' @return Data frame with information on allowed binding regions.
#' @keywords internal
retrieve.leader.region <- function(template.df, direction = c("fw", "rev"), start, end, gap.char) {
    if (length(start) == 0 || length(end) == 0) {
        return(NULL)
    }
    if (length(direction) == 0) {
        stop("Please provide the 'direction' argument.")
    }
    direction <- match.arg(direction)
    if (length(start) == 1 && length(end) == 1) {
        # start and end refer to intervals (uniform range)
        start <- rep(start, nrow(template.df))
        end <- rep(end, nrow(template.df))
    } # otherwise: regions are already ok in the input.
    if (direction == "rev") {
        # get reverse start and ends in absolute values:
        ori.start <- start
        start <- nchar(template.df$InputSequence) - end + 1
        end <- nchar(template.df$InputSequence) - ori.start + 1
    }
    # store start/end in possibly gappy input
    ali.start <- start
    ali.end <- end
    # gap count from start of seq until the end of the allowed region
    # adjust binding region in the ungapped sequence according to the number of gaps in between
    new.interval <- binding.interval.ungapped(template.df, start, end, gap.char)
    # work with start/end in ungapped seqs
    start <- new.interval[[1]]
    end <- new.interval[[2]]
    # use the ungapped sequences now to work on the correct positions 
    leader <- substring(template.df$Sequence, start, end)
    # deal with incorrect positions
    idx <- which(start < 0 & end < 0)  # interval is not in range -> no leader available
    if (length(idx) != 0) {
        start[idx] <- NA
        end[idx] <- NA
    }
    idx <- which(end < 0)  # only end is too small
    if (length(idx) != 0) {
        end[idx] <- 0  # no leader available
    }
    idx <- which(start < 0)  # only start is too early
    if (length(idx) != 0) {
        start[idx] <- 1
    }
    idx <- which(unlist(lapply(seq_along(end), function(x) end[x] > (nchar(template.df$Sequence[x]) + 
        1) && start[x] > (nchar(template.df$Sequence[x]) + 1))))
    if (length(idx) != 0) {
        # not in the defined range
        end[idx] <- NA
        start[idx] <- NA
    }
    idx <- which(unlist(lapply(seq_along(start), function(x) start[x] > (1 + nchar(template.df$Sequence[x])))))  #only start is larger than the range
    if (length(idx) != 0) {
        start[idx] <- 0
        end[idx] <- 0
    }
    idx <- which(unlist(lapply(seq_along(end), function(x) end[x] > (1 + nchar(template.df$Sequence[x])))))  #only end is larger than the range
    if (length(idx) != 0) {
        end[idx] <- nchar(template.df$Sequence)[idx]
    }
    df <- build_leader_df(direction, leader, start, end, ali.start, ali.end)
    return(df)
}
#' Building of Leader Data Frame.
#'
#' Constructs the leader data frame.
#'
#' @param direction The primer direction for which we are annotating binding regions.
#' @param leader The binding region sequence.
#' @param start The start positions of the binding region.
#' @param end The end positions of the binding region.
#' @param ali.start The start positions of the binding region in the aligned input.
#' @param ali.end The end positions of the binding region in the aligned input.
#' @return A data frame with binding region information.
#' @keywords internal
build_leader_df <- function(direction = c("fw", "rev"), leader, start, end,
                            ali.start, ali.end) {
    if (length(direction) == 0) {
        stop("Please provide the 'direction' argument.")
    }
    direction <- match.arg(direction)
    # initial: stored such that binding range can be modified later in relation to initial configuration
    df <- data.frame(Allowed_fw = character(length(leader)),  # binding region string
                    Allowed_rev = character(length(leader)),
                    # binding region intervals:
                    Allowed_Start_fw = numeric(length(start)), 
                    Allowed_Start_rev = numeric(length(start)),
                    Allowed_End_fw = numeric(length(end)), 
                    Allowed_End_rev = numeric(length(end)),
                    # the binding region in gapped sequence positions from alignemnts
                    Allowed_Start_fw_ali = numeric(length(start)), 
                    Allowed_Start_rev_ali = numeric(length(start)),
                    Allowed_End_fw_ali = numeric(length(end)), 
                    Allowed_End_rev_ali = numeric(length(end)),
                    # initial: the initially specified binding region (without gaps)
                    Allowed_Start_fw_initial = numeric(length(start)), 
                    Allowed_Start_rev_initial = numeric(length(start)), 
                    Allowed_End_fw_initial = numeric(length(end)), 
                    Allowed_End_rev_initial = numeric(length(end)),
                    # initial positions in aligned input:
                    Allowed_Start_fw_initial_ali = numeric(length(start)), 
                    Allowed_Start_rev_initial_ali = numeric(length(start)), 
                    Allowed_End_fw_initial_ali = numeric(length(end)), 
                    Allowed_End_rev_initial_ali = numeric(length(end)),
                    stringsAsFactors = FALSE)
    if (direction == "rev") {
        df$Allowed_rev <- leader
        df$Allowed_Start_rev <- start
        df$Allowed_End_rev <- end
        df$Allowed_Start_rev_ali <- ali.start
        df$Allowed_End_rev_ali <- ali.end
        df$Allowed_Start_rev_initial <- df$Allowed_Start_rev
        df$Allowed_End_rev_initial <- df$Allowed_End_rev
        df$Allowed_Start_rev_initial_ali <- df$Allowed_Start_rev_ali
        df$Allowed_End_rev_initial_ali <- df$Allowed_End_rev_ali
    } else {
        df$Allowed_fw <- leader
        df$Allowed_Start_fw <- start
        df$Allowed_End_fw <- end
        df$Allowed_Start_fw_ali <- ali.start
        df$Allowed_End_fw_ali <- ali.end
        df$Allowed_Start_fw_initial <- df$Allowed_Start_fw
        df$Allowed_End_fw_initial <- df$Allowed_End_fw
        df$Allowed_Start_fw_initial_ali <- df$Allowed_Start_fw_ali
        df$Allowed_End_fw_initial_ali <- df$Allowed_End_fw_ali
    }
    # we return only the modified binding direction such that 
    # previous settings in template.df are not overwritten:
    sel.cols <- grepl(direction, colnames(df))
    df <- df[,sel.cols]
    return(df)
}
#' Uniform Binding Ranges.
#'
#' Creates uniform binding regions for all templates.
#'
#' @param fw.interval Binding region for forward templates.
#' @param rev.interval Binding region for reverse templates.
#' @param template.df Template data frame.
#' @param gap.char The character for gaps in alignments.
#' @return Data frame with binding region information.
#' @keywords internal
create.uniform.leaders <- function(fw.interval, rev.interval, template.df, gap.char) {
    if (length(fw.interval) == 0 && length(rev.interval) == 0) {
        warning("No fw/rev intervals were specified for specifying a uniform binding range.")
    }
    if (length(fw.interval) == 2) {
        if (fw.interval[1] > fw.interval[2]) {
            stop("Forward binding region was no proper interval.")
        }
    }
    if (length(rev.interval) == 2) {
        if (rev.interval[1] > rev.interval[2]) {
            stop("Reverse binding region was no proper interval.")
        }
    }
    fw.leaders <- retrieve.leader.region(template.df, "fw", fw.interval[1], fw.interval[2], gap.char)
    rev.leaders <- retrieve.leader.region(template.df, "rev", rev.interval[1], rev.interval[2], gap.char)
    if (length(fw.leaders) == 0) {
        result <- rev.leaders
    } else if (length(rev.leaders) == 0) {
        result <- fw.leaders
    } else {
        result <- cbind(fw.leaders[,grep("_fw", colnames(fw.leaders))], 
                        rev.leaders[, grep("_rev", colnames(rev.leaders))])
    }
    col.order <- c("Allowed_Start_fw", "Allowed_End_fw", "Allowed_Start_rev", "Allowed_End_rev", 
        "Allowed_fw", "Allowed_rev")
    all.cols <- c(col.order[col.order %in% colnames(result)], setdiff(colnames(result), 
        col.order))
    result <- result[, all.cols]
    return(result)
}
#' Add Uniform Binding Regions.
#'
#' Augments a template data frame with uniform binding regions.
#'
#' @param lex.seq Template data frame
#' @param leaders Data frame with uniform binding regions.
#' @return Template data frame with updated binding regions.
#' @keywords internal
add.uniform.leaders.to.seqs <- function(lex.seq, leaders) {
    result <- lex.seq
    for (i in colnames(leaders)) {
        result[, i] <- leaders[, i]
    }
    return(result)
}
#' Read Individual Binding Regions
#'
#' Reads individual binding regions into a data frame.
#'
#' @param fasta.file Path to a FASTA file with binding regions.
#' @param direction String identifying whether the FASTA file
#' contains information pertaining to the binding region of forward
#' or reverse primers.
#' @param rm.keywords A vector of keywords that are used to remove templates whose headers contain any of the keywords.
#' @param gap.character The character for indicating gaps in sequences.
#' @return A data frame with individual binding regions.
#' @keywords internal
read.leaders <- function(fasta.file, direction = c("fw", "rev"), rm.keywords = NULL, gap.character) {
    if (length(fasta.file) == 0) {
        return(NULL)
    }
    d <- read.sequences(fasta.file, gap.character)
    if (direction == "fw") {
        colnames(d)[which(colnames(d) == "Sequence")] <- "Allowed_fw"
    } else {
        colnames(d)[which(colnames(d) == "Sequence")] <- "Allowed_rev"
    }
    if (length(rm.keywords) != 0) {
        d <- remove.seqs.by.keyword(d, rm.keywords)
    }
    return(d)
}
#' Parser for IMGT Groups.
#'
#' Parses IMGT group information contained in FASTA headers.
#'
#' @param IDs Group information strings to be parsed.
#' @return Data frame with structured group information.
#' @keywords internal
parse.IMGT.gene.groups <- function(IDs) {
    # regexp: (gene)-(subgroup)-(localization in the locus)-(extracted_labels)
    exp <- "([^-]+)-(.+)(\\*[^\\|]+)"
    # some genes are ORs (orphons) -> still the same gene group?  orphons appear
    # outside the main genetic loci but are still assigned to groups ..  clan
    # annotations: clan1: ighv1, ighv5, ighv7
    idx <- stringr::str_match(IDs, exp)
    if (any(apply(idx, 1, is.na))) {
        # at least one entry not in IMGT format
        return(NULL)
    }
    # remove OR annotation from GeneGroup
    colnames(idx) <- c("Full", "Gene", "Subgroup", "Allele")
    s <- strsplit(idx[, "Gene"], split = "/")
    g <- sapply(s, function(x) x[1])
    idx[, "Gene"] <- g
    # remove localization from Subgroup
    s <- strsplit(idx[, "Subgroup"], "Subgroup", split = "-")
    g <- sapply(s, function(x) x[1])
    idx[, "Subgroup"] <- g
    identifier <- paste(idx[, "Gene"], idx[, "Subgroup"], sep = "-")
    result <- data.frame(Group_Identifier = identifier, Gene_Group = idx[, "Gene"], 
        Sub_Group = idx[, "Subgroup"], Allele = idx[, "Allele"], 
        stringsAsFactors = FALSE)
    # message(result)
    return(result)
}
#' Unification of Leaders
#'
#' Unifies individual binding regions for forward and reverse primers.
#'
#' @param l.seq.fw Data frame with binding information for forward primers.
#' @param l.seq.rev Data frame with binding information for reverse primers.
#' @param lex.seq Template data frame.
#' @param gap.char The character for indicating alignment gaps.
#' @return Template data frame with annotated binding regions.
#' @keywords internal
unify.leaders <- function(l.seq.fw, l.seq.rev, lex.seq, gap.char) {
    fw.df <- NULL
    rev.df <- NULL
    if (length(l.seq.fw) == 0 && length(l.seq.rev) == 0) {
        return(lex.seq)
    }
    template.df <- lex.seq
    if (length(l.seq.fw) != 0) {
        fw.df <- get.leader.exon.regions.single(l.seq.fw, lex.seq, "fw", gap.char) 
        template.df <- add.uniform.leaders.to.seqs(template.df, fw.df)
    }
    if (length(l.seq.rev) != 0) {
        rev.df <- get.leader.exon.regions.single(l.seq.rev, lex.seq, "rev", gap.char)
        template.df <- add.uniform.leaders.to.seqs(template.df, rev.df)
    }
    template.df$Allowed_Start_fw_initial <- template.df$Allowed_Start_fw
    template.df$Allowed_End_fw_initial <- template.df$Allowed_End_fw
    template.df$Allowed_End_rev_initial <- template.df$Allowed_End_rev
    template.df$Allowed_Start_rev_initial <- template.df$Allowed_Start_rev
    template.df$Sequence <- ungap_sequence(template.df$InputSequence, gap.char) # use the ungapped seq from now on only
    return(template.df)
}
#' Assign Binding Regions
#' 
#' Augments a template data frame with individual binding regions.
#' 
#' @param lex.seqs Data frame with template sequences.
#' @param uni.leaders Data frame with individual allowed binding regions.
#' @return Template data frame with annotated binding regions.
#' @keywords internal
get.leader.exon.regions <- function(lex.seqs, uni.leaders) {
    if (length(uni.leaders) == 0) {
        return(NULL)
    }
    # add missing columns
    if (!"Allowed_Start_fw" %in% colnames(uni.leaders)) {
        uni.leaders$Allowed_Start_fw <- lex.seqs$Allowed_Start_fw
        uni.leaders$Allowed_End_fw <- lex.seqs$Allowed_End_fw
        uni.leaders$Allowed_fw <- ""
        uni.leaders$Allowed_Start_fw_initial <- uni.leaders$Allowed_Start_fw
        uni.leaders$Allowed_End_fw_initial <- uni.leaders$Allowed_End_fw
        
    }
    if (!"Allowed_Start_rev" %in% colnames(uni.leaders)) {
        uni.leaders$Allowed_Start_rev <- lex.seqs$Allowed_Start_rev
        uni.leaders$Allowed_End_rev <- lex.seqs$Allowed_End_rev
        uni.leaders$Allowed_rev <- ""
        uni.leaders$Allowed_Start_rev_initial <- uni.leaders$Allowed_Start_rev
        uni.leaders$Allowed_End_rev_initial <- uni.leaders$Allowed_End_rev
    }
    # add entries from lex.seqs, where uni.leaders didn't have anything
    fw.idx <- which(uni.leaders$Allowed_fw == "")
    rev.idx <- which(uni.leaders$Allowed_rev == "")
    if (length(fw.idx) != 0) {
        uni.leaders$Allowed_Start_fw[fw.idx] <- lex.seqs$Allowed_Start_fw[fw.idx]
        uni.leaders$Allowed_End_fw[fw.idx] <- lex.seqs$Allowed_End_fw[fw.idx]
        uni.leaders$Allowed_Start_fw_initial[fw.idx] <- lex.seqs$Allowed_Start_fw[fw.idx]
        uni.leaders$Allowed_End_fw_initial[fw.idx] <- lex.seqs$Allowed_End_fw[fw.idx]
    }
    if (length(rev.idx) != 0) {
        uni.leaders$Allowed_Start_rev[rev.idx] <- lex.seqs$Allowed_Start_rev[rev.idx]
        uni.leaders$Allowed_End_rev[rev.idx] <- lex.seqs$Allowed_End_rev[rev.idx]
        uni.leaders$Allowed_Start_rev_initial[rev.idx] <- lex.seqs$Allowed_Start_rev[rev.idx]
        uni.leaders$Allowed_End_rev_initial[rev.idx] <- lex.seqs$Allowed_End_rev[rev.idx]
        
    }
    # remove duplicate cols
    dup.idx <- which(duplicated(colnames(uni.leaders)))
    if (length(dup.idx) != 0) {
        uni.leaders <- uni.leaders[, -dup.idx]
    }
    return(uni.leaders)
}
#' Individual Binding Annotation
#'
#' Annotate individual binding regions.
#'
#' @param l.seq Data frame with individual binding regions.
#' @param lex.seq Template data frame.
#' @param direction The primer direction for which the binding info is valid.
#' @param gap.char The character for gaps in alignments.
#' @return Template data  frame with annotated binding regions.
#' @keywords internal
get.leader.exon.regions.single <- function(l.seq, lex.seq, 
                     direction = c("fw", "rev"), gap.char) {
    if (length(direction) == 0) {
        stop("Please provide the 'direction' argument.")
    }
    direction <- match.arg(direction)
    # match leader sequences to lex seqs
    # substring here is a necessary hack for mapping IMGT leaders to seqs
    idx <- lapply(lex.seq$Header, function(x) grep(substring(x, 1, 30), l.seq$Header, fixed = TRUE))  # hit from exon seqs to leaders
    # ensure that each template has a unique leader
    hit.len <- sapply(idx, length)
    nbr.of.hits <- sum(hit.len)
    if (any(hit.len == 0 | hit.len > 1)) {
        if (any(hit.len == 0)) {
            error.msg <- "Some templates did not have a corresponding binding region. The binding regions of templates for which no allowed region was specified were not changed."
            my.idx <- which(hit.len == 0)
            error.msg <- paste(error.msg, "The IDs of the sequences are:\n", paste(lex.seq$ID[my.idx], 
                collapse = "\n", sep = ""))
            my.warning("MissingLeaders", error.msg)
        } else {
            error.msg <- "Some templates had more than one corresponding binding region. Only the first specified binding region was considered in these cases."
            my.idx <- which(hit.len > 1)
            error.msg <- paste(error.msg, "The IDs of the sequences are:\n", paste(lex.seq$ID[my.idx], 
                collapse = "\n", sep = ""))
            my.warning("RedundantLeaders", error.msg)
        }
    } else if (nbr.of.hits != nrow(l.seq)) {
        # not all leaders were matched to a template sequence
        used.leaders <- unique(unlist(idx))
        not.used.leaders <- setdiff(seq_len(nrow(l.seq)), used.leaders)
        leader.info <- paste("Leaders with the following headers were not used: ", 
            paste(l.seq$Header[not.used.leaders], collapse = ", "), sep = "")
        my.warning("Not_all_leaders_matched", "Some allowed regions did not match any template header and were ignored.")
    }
    # select indices that are available if multiple leader regions are given for one
    # sample, take the 1st one arbtirarily
    sel.idx.zero <- which(sapply(idx, length) == 0)
    idx <- sapply(seq_along(idx), function(x) if (x %in% sel.idx.zero) {
        NA
    } else {
        idx[[x]][1]
    })
    sel.idx.lex <- which(!is.na(idx))  # ids of the templates where a leader seq was found
    idx <- unlist(idx)
    if (length(sel.idx.lex) == 0) {
        my.error("Leaders_no_matches", "No matches between templates and allowed binding regions. No binding regions were annotated")
    }
    sel.idx <- idx[sel.idx.lex]  # ids of the leader.seqs with a template
    ###### 
    template.df <- lex.seq[sel.idx.lex, ]
    idx.add <- which(!(colnames(l.seq) %in% colnames(template.df)))
    leader.col <- paste("Allowed", "_", direction, sep = "")
    template.df[, leader.col] <- l.seq$Allowed[sel.idx]  # replace the leader
    # rename the other columns that are already present and add
    cols <- setdiff(colnames(l.seq), leader.col)
    idx.replace <- which(cols %in% colnames(template.df))
    if (length(idx.replace) != 0) {
        template.df[, cols[idx.replace]] <- l.seq[sel.idx, cols[idx.replace]]
    }
    # add the new columns also ..
    if (length(idx.add) != 0) {
        template.df[, colnames(l.seq)[idx.add]] <- l.seq[sel.idx, idx.add]
    }
    
    # find leader position in sequences with exact matching if we have a forward
    # restriction, select the first match we find, otherwise select the last match we
    # find
    leader.pos <- lapply(seq_len(nrow(template.df)), function(x) gregexpr(template.df[x, leader.col], 
        template.df$InputSequence[x]))
    leader.idx <- NULL
    if (direction == "fw") {
        # always the take match we find
        leader.idx <- rep(1, length(leader.pos))
    } else {
        # always take the last match we find
        leader.idx <- unlist(lapply(seq_along(leader.pos), function(x) length(unlist(leader.pos[[x]]))))
    }
    leader.s <- sapply(seq_along(leader.pos), function(x) as.numeric(unlist(leader.pos[[x]])[leader.idx[x]]))
    leader.e <- sapply(seq_along(leader.pos), function(x) leader.s[x] + attr(leader.pos[[x]][[1]], 
        "match.length")[leader.idx[x]] - 1)
    # verify that all leader.s and leader.e are defined
    leader.s.ok <- leader.s > 0
    leader.e.ok <- leader.e > 0
    if (any(!leader.s.ok) || any(!leader.e.ok)) {
        warn.msg <- "Some of the input allowed binding sites could not be found in the templates. The affected sequences are:\n"
        idx1 <- which(!leader.s.ok)
        idx2 <- which(!leader.e.ok)
        idx <- union(idx1, idx2)
        warn.msg <- paste(warn.msg, paste(template.df$ID[idx], collapse = "\n", sep = ""))
        my.warning("LeadersNotFound", warn.msg)
    }
    final.df <- retrieve.leader.region(template.df, direction, leader.s, leader.e, gap.char) # retrieve leader region: overwrites values of other directions
    if (FALSE) {
        if (any(leader.s.ok)) {
        # set allowed region start/end
        allowed.col.e <- paste("Allowed_End", "_", direction, sep = "")
        allowed.col.s <- paste("Allowed_Start", "_", direction, sep = "")
        template.df[leader.s.ok, allowed.col.s] <- leader.s[leader.s.ok]
        template.df[leader.s.ok, allowed.col.e] <- leader.e[leader.s.ok]
        }
    }
    return(final.df)
}
#' Storing Templates to Disk.
#'
#' Stores a set of templates as a FASTA or CSV file.
#'
#' @param template.df An object of class \code{Templates} to be stored to disk.
#' @param fname The path to the file where the templates should be stored.
#' @param ftype A character vector giving the filetype.
#' This can either be "FASTA" or "CSV" (default: "FASTA").
#' @return Stores templates to \code{fname}.
#' @family template functions
#' @keywords Templates
#' @export
#' @examples
#' data(Ippolito)
#' # Store templates as FASTA
#' fname.fasta <- tempfile("my_templates", fileext = ".fasta")
#' write_templates(template.df, fname.fasta)
#' # Store templates as CSV
#' fname.csv <- tempfile("my_templates", fileext = ".csv")
#' write_templates(template.df, fname.csv, "CSV")
write_templates <- function(template.df, fname, ftype = c("FASTA", "CSV")) {
    ftype <- match.arg(ftype)
    if (ftype == "FASTA") {
        seqs <- template.df$Sequence
        labels <- template.df$ID
        seqinr::write.fasta(as.list(seqs), as.list(labels), fname)
    } else if (ftype == "CSV") {
        write.csv(template.df, file = fname, row.names = FALSE)
    } else {
        warning("Unknown filetype: ", ftype)
    }
}

#' Scoring of Template Conservation.
#'
#' Determines the sequence conservation scores of a set of templates
#' using Shannon entropy.
#'
#' @param template.df A \code{Templates} object providing the sequence
#' conservation shall be determined.
#' @param gap.char The alignment gap character. By default, this is set to "-".
#' @param win.len The size of a window for evaluating conservation.
#' The default window size is set to 30.
#' @param by.group Whether the determination of binding regions 
#' should be stratified according to the groups defined in \code{template.df}.
#' The default is \code{TRUE}.
#' @return A list containing \code{Entropies} and \code{Alignments}.
#' Entropies are given as a data frame with conservation scores. 
#' Each column indicates a position in the alignment of template sequences
#' and each row gives the entropies of the sequences 
#' belonging to a specific group of template sequences.
#' Alignments are given as lists of \code{DNABin} objects, where each
#' object gives the alignment corresponding to one group of template sequences.
#' @export
#' @note Requires the MAFFT software for multiple alignments (http://mafft.cbrc.jp/alignment/software/).
#' @examples
#' \dontrun{
#' data(Ippolito)
#' entropy.data <- score_conservation(template.df, gap.char = "-", win.len = 18, by.group = TRUE)
#' }
score_conservation <- function(template.df, gap.char = "-", 
                               win.len = 30, by.group = TRUE) {
    
    # only align sequences if not aligned already
    if (!check.tool.function()["MAFFT"]) {
        stop("MAFFT is required for scoring the conservation (http://mafft.cbrc.jp/alignment/software/).")
    }
    if (!any(grepl(gap.char, template.df$InputSequence))) {
        # no gap char found -> sequences are not aligned yet
        message("Aligning sequences ...")
        ali <- align.seqs(template.df$InputSequence, template.df$ID)
    } else {
        message("Assuming that the sequences are pre-aligned.")
        ali <- seqinr::as.alignment(nb = length(template.df$InputSequence),
                                    nam = template.df$ID,
                                    seq = template.df$InputSequence)
    }
    if (is.null(ali)) {
        # alignment wasn't possible
        return(NULL)
    }
    # sanitize the sequences: (carriage returns cause problems in windows)
	ali$seq <- gsub("\r", "", ali$seq)
    bin <- ape::as.DNAbin(ali)
    # create groups in the alignment according to user-anntotation
    # may automatic grouping in the future (clustering?)
    if (by.group) {
        groups <- unique(template.df$Group)
        bins <- lapply(seq_along(groups), function(x) bin[which(template.df$Group == groups[x]),])
    } else {
        groups <- "default"
        bins <- list(bin)
    }
    # determine conservation according to Shannon entropy for each alignment group
    ali.entropy <- lapply(bins, function(x) shannon.entropy(as.character(x)))
    # create a data frame representation
    entropy.df <- do.call(rbind, ali.entropy)
    rownames(entropy.df) <- groups
    colnames(entropy.df) <- seq_len(ncol(entropy.df))
    names(bins) <- groups
    out <- list("Alignments" = bins, "Entropies" = entropy.df)
    return(out)
}
#' Identification of Gappy Columns in Alignments.
#' @param bins A list of \code{DNABin} alignments.
#' @param gap.cutoff The required percentage of gaps
#' for consideration as a gap column.
#' @param gap.char The gap character in the alignments.
#' @return A list with indices giving the gap columns
#' for every alignment in \code{bins}.
#' @keywords internal
detect.gap.columns <- function(bins, gap.cutoff = 0.95, gap.char = "-") {
    gap.cutoff <- 0.95 # if there is a column 95% gaps, ignore the whole region
    gap.freq <- lapply(bins, function(x) {
        mat <- as.character(x)
        count <- unlist(lapply(seq_len(ncol(mat)), function(y) 
            length(which(mat[,y] == gap.char))))
        freq <- count / nrow(mat)
    })
    # exclude regions with gap percentage above cutoff -> useless to construct primers here 
    mask.idx <- lapply(gap.freq, function(x) which(x >= gap.cutoff))
    return(mask.idx)
}
#' Plot of Template Sequence Conservation.
#'
#' Plots the template sequence conservation (range [0,1]) according to
#' the Shannon entropy of the sequences. 
#'
#' @param entropy.df A data frame with entropies.
#' Each row gives the entropies of a group of related template
#' sequences for all columns of the alignment.
#' @param alignments A list with \code{DNABin} alignment objects
#' corresponding to the groups (rows) in the alignment.
#' @param template.df The \code{Templates} object for which 
#' the conservation has been determined.
#' @param gap.char The gap char in the alignments. By default,
#' \code{gap.char} is set to "-".
#' @return A plot showing the degree of sequence conservation in the templates.
#' @export
#' @note Computing the conservation scores for the plot requires the MAFFT software for multiple alignments (http://mafft.cbrc.jp/alignment/software/).
#' @examples
#' \dontrun{
#' data(Ippolito)
#' # Select binding regions for every group of templates and plot:
#' template.df <- select_regions_by_conservation(template.df, win.len = 30)
#' plot_conservation(attr(template.df, "entropies"), attr(template.df, "alignments"), template.df)
#' # Select binding regions for all templates and plot:
#' data(Ippolito)
#' template.df <- select_regions_by_conservation(template.df, by.group = FALSE)
#' plot_conservation(attr(template.df, "entropies"), attr(template.df, "alignments"), template.df)
#' } 
plot_conservation <- function(entropy.df, alignments, template.df, gap.char = "-") {
    # set gappy columns to 0 conservation for all groups
    gap.idx <- detect.gap.columns(alignments, gap.cutoff = 0.95, gap.char = gap.char)
    m <- match(rownames(entropy.df), names(gap.idx))
    for (i in seq_len(nrow(entropy.df))) { # for each alignment group
        entropy.df[i, gap.idx[[m[i]]]] <- 1 # max entropy
    }
    plot.df <- reshape2::melt(entropy.df, value.name = "Entropy")
    # since we have [0,1] interval entropies, we can simply convert to conservation:
    colnames(plot.df) <- c("Group", "Position", "Entropy")
    plot.df$Conservation <- 1 - plot.df$Entropy
    # get region data for every group
    get.region.df <- function(allowed.start.fw.old, allowed.end.fw.old,
                              allowed.start.fw.new, allowed.end.fw.new,
                              allowed.start.rev.old, allowed.end.rev.old,
                              allowed.start.rev.new, allowed.end.rev.new) {
        # select unique interval here, but how?
        fw.ori.range <- unique(cbind(allowed.start.fw.old,
                                allowed.end.fw.old))
        fw.ori.df <- data.frame("Start" = fw.ori.range[,1],
                                "End" = fw.ori.range[,2],
                                "Direction" = "fw",
                                "Type" = "Old")
        fw.new.range <- unique(cbind(allowed.start.fw.new,
                                allowed.end.fw.new))
        fw.new.df <- data.frame("Start" = fw.new.range[,1],
                                "End" = fw.new.range[,2],
                                "Direction" = "fw",
                                "Type" = "New")
        rev.ori.range <- unique(cbind(allowed.start.rev.old,
                                allowed.end.rev.old))
        rev.ori.df <- data.frame("Start" = rev.ori.range[,1],
                                "End" = rev.ori.range[,2],
                                "Direction" = "rev",
                                "Type" = "Old")
        rev.new.range <- unique(cbind(allowed.start.rev.new,
                                allowed.end.rev.new))
        rev.new.df <- data.frame("Start" = rev.new.range[,1],
                                "End" = rev.new.range[,2],
                                "Direction" = "rev",
                                "Type" = "New")
        range.df <- rbind(fw.ori.df, fw.new.df, 
                          rev.ori.df, rev.new.df)
        return(range.df)
    }
    groups <- names(alignments)
    allowed.df <- vector("list", length(groups))
    for (i in seq_along(alignments)) {
        if (groups[i] == "default") { 
            # capture all templates
            sub.df <- template.df
        } else {
            # select individual groups
            sub.df <- template.df[template.df$Group == groups[i],]
        }
        region.df <- get.region.df(sub.df$Allowed_Start_fw_initial_ali,
                               sub.df$Allowed_End_fw_initial_ali,
                               sub.df$Allowed_Start_fw_ali,
                               sub.df$Allowed_End_fw_ali,
                               sub.df$Allowed_Start_rev_initial_ali,
                               sub.df$Allowed_End_rev_initial_ali,
                               sub.df$Allowed_Start_rev_ali,
                               sub.df$Allowed_End_rev_ali)
        region.df$Group <- groups[i]
        allowed.df[[i]] <- region.df
    }
    allowed.df <- do.call(rbind, allowed.df)
    group.colors <- RColorBrewer::brewer.pal(8, "Set1")
    group.colors <- grDevices::colorRampPalette(group.colors)(length(groups))
    group.colors <- c(group.colors, "grey40")
    # only plot the new allowed regions
    allowed.new <- allowed.df[allowed.df$Type == "New",]
    allowed.old <- allowed.df[allowed.df$Type == "Old",]
    # get max bar height for any group for a position
    #max.conservation <- plyr::ddply(plot.df, c("Position"), plyr::summarize, "TotalConservation" = sum(substitute(Conservation)))
    #box.extent <- max(max.conservation$TotalConservation) * 0.1
    box.extent <- max(plot.df$Conservation) * 0.1
    box.ymax <- -0.05 * max(plot.df$Conservation)
    box.ymin <- box.ymax - box.extent
    # select only unique old regions
    allowed.old <- allowed.old[row.names(unique(allowed.old[, c("Start", "End")])), ]
    p <- ggplot() + 
        geom_bar(data = plot.df, aes_string(x = "Position", 
                 y = "Conservation",
                 fill = "Conservation"),
                 stat = "identity",
                 position = "stack") +
        xlab("Position") +
        ylab("Conservation") +
        scale_y_continuous(labels = scales::percent) +
        facet_wrap(~Group, ncol = 2) 
    # add rectangles for binding region annotation
    p <- p + geom_rect(data = allowed.new, 
                      aes_string(xmin = "Start", 
                                 xmax = "End", 
                                 ymax = box.ymax,
                                 ymin = box.ymin),
                                 #color = "Type",
                                 alpha = 1) +
            # don't plot old regions: may not be updated to aligned binding regions
            #geom_rect(data = allowed.old,
                #aes_string(xmin = "Start", 
                                 #xmax = "End", 
                                 #fill = "Type",
                                 #ymax = box.ymax,
                                 #ymin = box.ymin),
                                 #alpha = 0.2) +
      #scale_fill_manual(values = group.colors) +
      scale_fill_gradient()
      ggtitle("Sequence conservation")
    return(p)
}
#' Updates Binding Region in the Alignment by conservation.
#' @param template.df A \code{Templates} object.
#' @param bins A list with \code{DNAbin} alignments, one for each group of template sequences.
#' @param entropy.df A data frame with entropy information.
#' @param gap.char The gap character for alignments.
#' @param win.len The desired length of the new binding region.
#' @param direction The direction for which the binding range shall be adjusted.
#' @return A \code{Templates} object with modified binding regions.
#' @keywords internal
update.binding.ranges.by.conservation <- function(template.df, 
                                    bins,
                                    entropy.df,
                                    gap.char = "-", 
                                    win.len = 30, direction = c("fw", "rev")) {
    if (length(direction) == 0) {
        stop("Please provide the 'direction' argument.")
    }
    direction <- match.arg(direction)
    if (direction == "fw") {
        allowed.region <- cbind(template.df$Allowed_Start_fw_initial_ali, 
                                template.df$Allowed_End_fw_initial_ali)
    } else {
        allowed.region <- cbind(template.df$Allowed_Start_rev_initial_ali, 
                                template.df$Allowed_End_rev_initial_ali)

    }
    step.size <- 1 # shift each window by 5 positions for speedup of computations
    primer.range <- create.primer.ranges(allowed.region[,2], rep(win.len, nrow(allowed.region)), 
                                        start.position = allowed.region[,1], step.size = step.size,
                                        groups = template.df$Group) 
    primer.ranges <- parallel::mclapply(seq_along(bins), 
        function(x) select.primer.region.by.conservation(primer.range[primer.range$Group == names(bins)[x], ], 
            entropy.df[x,], 0.1, bins[[x]]) # select top 10% within the allowed region
    )
    # annotate ranges with groups
    for (i in seq_along(primer.ranges)) {
        primer.ranges[[i]] <- cbind(Group = names(bins)[i], primer.ranges[[i]]) # TODO error
    }
    ranges <- do.call(rbind, primer.ranges)
    # select range with smallest entropy for every group
    selected.ranges <- plyr::ddply(ranges, "Group", function(x) plyr::arrange(x, substitute(Entropy))[1,])
    # assign new binding regions for every group of template sequences individually
    new.templates <- template.df
    for (i in seq_along(selected.ranges$Group)) {
        group <- selected.ranges$Group[i]
        ali <- as.character(bins[[i]])
        ali.seqs <- unlist(lapply(seq_len(nrow(ali)),function(x) paste0(ali[x,], collapse = "")))
        # identify templates belonging to this group
        idx <- which(template.df$Group == group)
        sub.df <- template.df[idx,]
        # set the new aligned sequence
        sub.df$InputSequence <- ali.seqs
        region <- c(selected.ranges[i, "Start"], selected.ranges[i, "End"])
        if (direction == "fw")  {
            region.fw <- region
            region.rev <- NULL
        } else {
            region.fw <- NULL
            # convert from absolute posi to relative posi from the end of the seq
            region.rev <- rev(ncol(ali) - region)
        }
        sub.df <- assign_binding_regions(sub.df, region.fw, region.rev)
        # store old (initial binding regions)
        if (direction == "fw") {
            sub.df$Allowed_Start_fw_initial_ali <- template.df[idx, "Allowed_Start_fw_initial_ali"]
            sub.df$Allowed_End_fw_initial_ali <- template.df[idx, "Allowed_End_fw_initial_ali"]
        } else {
            sub.df$Allowed_Start_rev_initial_ali <- template.df[idx, "Allowed_Start_rev_initial_ali"]
            sub.df$Allowed_End_rev_initial_ali <- template.df[idx, "Allowed_End_rev_initial_ali"]
        }
        # overwrite template entry
        new.templates[idx,] <- sub.df
    }
    return(new.templates)
}

#' Selection of Primer Binding Regions by Conservation.
#'
#' Computes Shannon entropy for putative binding regions
#' and determines the most conserved regions.
#'
#' @param template.df A \code{Templates} object containing
#' the template sequences for which the binding regions
#' shall be determined according to conservation.
#' @param gap.char The alignment gap character. This is "-" by default.
#' @param win.len The extent of the desired primer binding region.
#' This should be smaller than the \code{allowed.region}. The default is 30.
#' @param by.group Shall the determination of binding regions be stratified
#' according to the groups defined in \code{template.df}. By default,
#' this is set to \code{TRUE}.
#' @param direction Whether regions shall be selected for primers of
#' both directions ("both"), forward primers ("fw"), or reverse primers ("rev").
#' The default is "both".
#' @return A \code{Templates} object with adjusted binding regions.
#' The attribute \code{entropies} gives a data frame with positional entropies
#' and the attribute \code{alignments} gives the alignments of the templates.
#' @export
#' @note Requires the MAFFT software for multiple alignments (http://mafft.cbrc.jp/alignment/software/).
#' @examples
#' \dontrun{
#' data(Ippolito)
#' new.template.df <- select_regions_by_conservation(template.df)
#' }
select_regions_by_conservation <- function(template.df, 
                                    gap.char = "-", 
                                    win.len = 30, by.group = TRUE,
                                    direction = c("both", "fw", "rev")) {
    if (!is(template.df, "Templates")) {
        stop("Please input an object of class 'Templates'.")
    }
    direction <- match.arg(direction)
    # TODO: save time by not computing the entropy of the full sequence, but only the entropy of the possible binding regions?
    # TODO: individual binding regions are not really respected as selection is via groups and not individual sequences.
    entropy.data <- score_conservation(template.df, gap.char, win.len, by.group)
    bins <- entropy.data$Alignments # list with DNABin objects
    entropy.df <- entropy.data$Entropies # list with data frames giving the entropies
    if (direction == "fw" || direction == "rev") {
        template.df <- update.binding.ranges.by.conservation(template.df, bins, entropy.df, gap.char, win.len, direction)
    } else if (direction == "both") {
        t.fw <- update.binding.ranges.by.conservation(template.df, bins, entropy.df, gap.char, win.len, "fw")
        template.df <- update.binding.ranges.by.conservation(t.fw, bins, entropy.df, gap.char, win.len, "rev")
    }
    # add additional data (for plotting) as attributes
    attr(template.df, "entropies") <- entropy.df
    attr(template.df, "alignments") <- bins
    return(template.df)
}

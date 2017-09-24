#' Output Functionalities.
#'
#' @rdname Output
#' @name Output
#' @description
#' \describe{
#' \item{\code{write_primers}}{Writes a set of primers to disk, either as a FASTA or CSV file.}
#' \item{\code{write_settings}}{ Stores primer analysis settings to a file in XML format.}
#' \item{\code{write_templates}}{Stores a set of templates as a FASTA or CSV file.}
#' \item{\code{create_report}}{Creates a PDF report for analyzed primer sets.}
#' \item{\code{create_coverage_xls}}{Creation of an XLS spreadsheet 
#' providing an overview of the covered
#' template sequences for each primer. Each cell in the spreadsheet
#' indicates a coverage event between a primer and template using
#' color codes. Identified coverage events are indicated by green, while
#' primer-template pairs without coverage are indicated by red.
#' In case that a primer binding condition (see \code{\link{CoverageConstraints}})
#' was active when computing the coverage, the numeric value of the
#' coverage condition is annotated for each cell.}
#' }
#' 
#' @param primer.df An object of class \code{Primers}.
#' @param fname The path to the output file.
#' @param ftype A character vector giving the type of the file.
#' This can either be "FASTA" or "CSV" (default: "FASTA").
#' @param template.df An object of class \code{Templates}.
#' @param settings A \code{DesignSettings} object to be stored to disk.
#' @param primers To create a report for a single primer set, please provide
#' an evaluated \code{Primers} object.
#' For creating a report comparing multiple primer sets, please provide
#' a list of \code{Primers} objects.
#' @param templates If \code{primers} is a \code{Primers} object, \code{templates} should be a \code{Templates} object.
#' If \code{primers} is a list of \code{Primers} objects, \code{templates}
#' should be a list of \code{Templates} objects of the same length as \code{primers}.
#' @param sample.name An identifier for your analysis. By default ( 
#' \code{NULL}), the sample identifier is selected from the
#' \code{Run} column of the input templates.
#' @param used.settings A named list (with fields \code{fw} and \code{rev}) containing the relaxed settings
#' for designing forward/reverse primers. By default (\code{NULL}), 
#' the relaxed settings are not shown in the report.
#' @param ... \code{required.cvg} (optional, default: 1), the desired coverage ratio if \code{primers} is a single primer set.
NULL
 

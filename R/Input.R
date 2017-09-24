#' Input Functionalities.
#'
#' @rdname Input
#' @name Input
#' @description
#' \describe{
#' \item{\code{read_primers}}{Reads one or multiple input files with 
#' primer sequences. The input can either be in FASTA or in CSV format.}
#' \item{\code{read_templates}}{Read one or multiple files with template sequences in FASTA or CSV format.}
#' \item{\code{read_settings}}{Loads primer analysis settings from an XML file.}
#' \item{\code{Templates}}{The \code{Templates} class encapsulates a data frame
#' containing the sequencs of the templates, their binding regions,
#' as well as additional information (e.g. template coverage).}
#' \item{\code{Primers}}{The \code{Primers} class encapsulates a data frame
#' representing a set of primers. Objects of this class
#' store all properties associated with a set of primers,
#' for example the results from evaluating the properties
#' of a primer set or from determining its coverage.}
#' }
#' @param fname Character vector providing either a single or multiple
#' paths to FASTA or CSV files.
#' @param fw.id For FASTA input, the identifier for forward primers in the FASTA headers.
#' @param rev.id For FASTA input, the identifier for reverse primers in the FASTA headers.
#' @param merge.ambig Indicates whether similar primers should be merged
#' ("merge") using IUPAC ambiguity codes or whether primers should 
#' be disambiguated ("unmerge").
#' By default \code{merge.ambig} is set to "none", leaving primers as they are.
#' @param max.degen A scalar numeric providing the maximum
#' allowed degeneracy for merging primers if
#' \code{merge.ambig} is set to "merge".
#' Degeneracy is defined by the number of disambiguated sequences 
#' that are represented by a degenerate primer.
#' @param template.df An object of class \code{Templates}.
#' If \code{template.df} is provided for \code{read_primers} then the primers are checked
#' for restriction sites upon input; otherwhise they are not checked.
#' @param adapter.action The action to be performed when \code{template.df} is
#' provided for identifying adapter sequences.
#' Either "warn" to issue warning about adapter sequences or
#' "rm" to remove identified adapter sequences. The default is "warn".
#' @param sample.name An identifier for the input primers. 
#' @param updateProgress A Shiny progress callback function. This is
#' \code{NULL} by default such that no progress is tracked.
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
#' @param gap.char The character in the input file representing gaps. 
#' Gaps are automatically removed upon input and the default character is "-".
#' @param run An identifier for the set of template sequences. By default,
#' \code{run} is \code{NULL} and its value is set via \code{template.file}.
#' @param filename Path to a valid XML file containing the 
#' primer analysis settings. By default, \code{filename} is set
#' to all settings that are shipped with openPrimeR and the lexicographically
#' first file is loaded.
#' @param frontend Indicates whether settings shall be loaded for the Shiny frontend. 
#' In this case no unit conversions for the PCR settings are performed.
#' The default setting is \code{FALSE} such that the correct units are used.
#' @param ... A data frame fulfilling the structural requirements
#' for initializing a \code{Templates} or \code{Primers} object.
NULL

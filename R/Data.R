#' Data Sets
#'
#' @rdname Data
#' @name Data
#' @docType data
#' @keywords datasets
#' @description
#' \describe{
#' \item{\code{Ippolito}}{IGHV primer data from Ippolito et al.}
#' \item{\code{Tiller}}{IGHV primer data from Tiller et al.}
#' \item{\code{Comparison}}{Evaluated primer sets targeting the 
#' functional human IGH immunoglobulin
#' genes. The sets were generated using the default evaluation settings
#' of openPrimeR. The primer sets were gathered from
#' IMGT and the literature.}
#' \item{\code{RefCoverage}}{Experimental results of multiplex PCR.}
#' }
#'  
#' @format For the \code{RefCoverage} data set, the \code{feature.matrix} data frame contains the properties
#' of the primer set from Tiller et al. as well as a primer set
#' that was designed by openPrimeR. The column \code{Experimental_Coverage}
#' indicates the experimentally determined coverage, while the other
#' columns relate to properties of the primers that were computed
#' with openPrimeR. 
#' The \code{ref.data} list contains the raw experimental coverage
#' of individual primers from the primer sets from Tiller and openPrimeR, which
#' both target templates from the IGH locus.
#' The rows of the data frames indicate primers and the columns indicate
#' IGH templates for which experimental coverage was determined.
#' The cell entries are hex codes. Each hex code represents a color 
#' indicating a certain experimental coverage status. 
#' Hex codes representing red shades indicate no or little amplification, 
#' while hex codes for green shades indicate high yields.
#'
#' For the \code{Ippolito} data set, 
#' \code{primer.df} provides a \code{Primers} object containing
#' the evaluated set of primers from Tiller et al.
#' \code{template.df} provides
#' a \code{Templates} object containing functional, human IGHV
#' templates for, and \code{settings} provides a 
#' \code{DesignSettings} object providing the used analysis settings.
#'
#' For the \code{Comparison} data set, \code{primer.data} and \code{template.data} are
#' lists of \code{Primers} and \code{Templates} objects, respectively.
#' 
#' For the \code{Tiller} data set, \code{tiller.primer.df} 
#' provides a \code{Primers} object, \code{tiller.template.df} provides 
#' the corresponding \code{Templates} object, and \code{tiller.settings} 
#' provides the \code{DesignSettings} object that was used
#' for evaluating \code{tiller.primer.df}.
NULL

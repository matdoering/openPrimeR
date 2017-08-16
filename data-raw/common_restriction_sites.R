###
# generate common restriction site Rdata
####
# seqRFLP data from REBASE:
# pkg: seqRFLP (>= 1.0.1)
utils::data("enzdata", package = "seqRFLP")
# common restriction enzymes retrieved from:
# https://www.addgene.org/mol-bio-reference/restriction-enzymes/
# annotate common sites in enzdata
common.file <- "src/analyses/common_restriction_sites.txt"
common.sites <- as.character(read.table(common.file)[,1])
m <- match(common.sites, enzdata$nam)
enzdata$common <- FALSE
enzdata$common[m[!is.na(m)]] <- TRUE
devtools::use_data(enzdata, pkg = "src/openPrimeR", internal = TRUE) # put this in the inst/extdata folder since it's for internal use by the pkg


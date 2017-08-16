# load templates
fasta.file <- system.file("extdata", "IMGT_data", "templates", 
              "Homo_sapiens_IGH_functional_exon.fasta", package = "openPrimeR")
hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")
template.df <- read_templates(fasta.file, hdr.structure, "|", "GROUP") 

# load primers
primer.location <- system.file("extdata", "IMGT_data", "primers", "IGHV", 
                               "Ippolito2012.fasta", package = "openPrimeR")
primer.df <- read_primers(primer.location, "_fw", "_rev")

# load settings
settings.xml <- system.file("extdata", "settings", 
                "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
settings <- read_settings(settings.xml)

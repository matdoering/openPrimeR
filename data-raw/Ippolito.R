# Store Ippolito's primers for IGH, the templates, and the settings.
fasta.file <- system.file("extdata", "IMGT_data", "templates", 
            "Homo_sapiens_IGH_functional_exon.fasta", 
            package = "openPrimeR")
hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")
# load templates
template.df <- openPrimeR::read_templates(fasta.file, hdr.structure, "|", "GROUP", rm.keywords = "partial")
# assign leader:
leader.file <- system.file("extdata", "IMGT_data", "templates",
            "Homo_sapiens_IGH_functional_leader.fasta", package = "openPrimeR")
template.df <- assign_binding_regions(template.df, leader.file, NULL)
primer.location <- system.file("extdata", "IMGT_data", "primers", "IGHV", 
                        "Ippolito2012.fasta", package = "openPrimeR")
primer.df <- openPrimeR::read_primers(primer.location, "_fw", "_rev")
# load settings for analysis
filename <- system.file("extdata", "settings", 
            "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
settings <- openPrimeR::read_settings(filename)
# evaluate constraints
# allow also other binding events:
conOptions(settings)$allowed_other_binding_ratio <- 1.0
conOptions(settings)$allowed_mismatches <- 5
primer.df <- openPrimeR::check_constraints(primer.df, template.df, settings)
template.df <- openPrimeR::update_template_cvg(template.df, primer.df)
out.loc <- system.file("data", "Ippolito.rda", package = "openPrimeR")
# compression level: 'xz' ensures small size
save(primer.df, template.df, settings, file = out.loc, compress = "xz")

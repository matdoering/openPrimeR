# Store Tiller's primers for IGH, the templates, and the settings.
fasta.file <- system.file("extdata", "IMGT_data", "templates", 
            "Homo_sapiens_IGH_functional_exon.fasta", 
            package = "openPrimeR")
hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")
# load templates
tiller.template.df <- openPrimeR::read_templates(fasta.file, hdr.structure, "|", "GROUP", rm.keywords = "partial")
# assign leader:
leader.file <- system.file("extdata", "IMGT_data", "templates",
            "Homo_sapiens_IGH_functional_leader.fasta", package = "openPrimeR")
tiller.template.df <- assign_binding_regions(tiller.template.df, leader.file, NULL)
# load primers
primer.location <- system.file("extdata", "IMGT_data", "primers", "IGHV", 
                        "Tiller2008_1st.fasta", package = "openPrimeR")
tiller.primer.df <- openPrimeR::read_primers(primer.location, "_fw", "_rev")
# load settings for analysis
filename <- system.file("extdata", "settings", 
            "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
tiller.settings <- openPrimeR::read_settings(filename)
# evaluate constraints
# allow also other binding events:
conOptions(tiller.settings)$allowed_other_binding_ratio <- 1.0
conOptions(tiller.settings)$allowed_mismatches <- 5
tiller.primer.df <- openPrimeR::check_constraints(tiller.primer.df, tiller.template.df, tiller.settings)
tiller.template.df <- openPrimeR::update_template_cvg(tiller.template.df, tiller.primer.df)
out.loc <- file.path(system.file("data",  package = "openPrimeR"), "Tiller.rda")
# compression level: 'xz' ensures small size
save(tiller.primer.df, tiller.template.df, tiller.settings, file = out.loc, compress = "xz")

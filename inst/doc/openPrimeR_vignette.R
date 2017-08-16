## ----vignette_options, echo = FALSE, message = FALSE, warning = FALSE----
knitr::opts_chunk$set(collapse = TRUE, comment = "#>")
library(openPrimeR)
ggplot2::theme_set(ggplot2::theme_grey(base_size = 12)) 

## ----check_dependencies, message = FALSE, warning = FALSE, eval = FALSE----
#  library(openPrimeR)

## ----loading_data_table, echo=FALSE, message=FALSE, warnings=FALSE, results='asis'----
tabl <- "
| Task             | Templates | Primers      | Input file format  |
|------------------|-----------|--------------|--------------------|
| Design primers   | &#10003;  |              | FASTA              |
| Analyze primers  | &#10003;  | &#10003;     | FASTA, CSV         |
| Compare primers  | &#10003;  | &#10003;     | FASTA, CSV         |
"
cat(tabl) # output the table in a format good for HTML/PDF/docx conversion

## ----load_templates, message = FALSE, warning = FALSE--------------------
fasta.file <- system.file("extdata", "IMGT_data", "templates", 
                "Homo_sapiens_IGH_functional_exon.fasta", package = "openPrimeR")
hdr.structure <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")
seq.df <- read_templates(fasta.file, hdr.structure, delim = "|", id.column = "GROUP")

## ----header_structure, message = FALSE, warning = FALSE------------------
seq.df$Header[1]

## ----simple_load_templates, message = FALSE, warning = FALSE-------------
seq.df.simple <- read_templates(fasta.file)

## ----assign_uniform_binding, message = FALSE, warning = FALSE------------
template.df.uni <- assign_binding_regions(seq.df, fw = c(1,50), rev = c(1,40))

## ----uniform_binding_regions, message = FALSE, warning = FALSE-----------
template.df.uni$Allowed_fw[1]
template.df.uni$Allowed_rev[1]

## ----assign_individual_binding, message = FALSE, warning = FALSE---------
l.fasta.file <- system.file("extdata", "IMGT_data", "templates", 
                "Homo_sapiens_IGH_functional_leader.fasta", package = "openPrimeR")
template.df <- assign_binding_regions(seq.df, fw = l.fasta.file, rev = NULL)

## ----assign_individual_binding_out, message = FALSE, warning = FALSE-----
template.df$Allowed_End_fw[1]
template.df$Allowed_End_fw[150]

## ----assign_individual_binding_example_rev, message = FALSE, warning = FALSE----
template.df$Allowed_rev[1]

## ----assign_optimized_binding, message = FALSE, warning = FALSE, eval = FALSE----
#  # requires ViennaRNA
#  template.df.uni.opti <- assign_binding_regions(seq.df, c(1,50), c(1,40),
#                              optimize.region = TRUE, primer.length = 20)

## ----assign_optimized_binding_adjusted, message = FALSE, warning = FALSE, eval = FALSE----
#  template.df.uni$Allowed_fw[1]
#  template.df.uni.opti$Allowed_fw[1]

## ----available_settings, message = FALSE, warning = FALSE----------------
list.files(system.file("extdata", "settings", package = "openPrimeR"), pattern = "*\\.xml")

## ----load_settings, message = FALSE, warning = FALSE---------------------
settings.xml <- system.file("extdata", "settings", 
                    "C_Taq_PCR_high_stringency.xml", package = "openPrimeR")
settings <- read_settings(settings.xml)

## ----cvg_constraint_setup, message = FALSE, warning = FALSE--------------
cvg_constraints(settings) <- list("primer_efficiency" = c("min" = 0.001))

## ----change_settings1, message = FALSE, warning = FALSE------------------
design.settings <- settings
constraints(design.settings) <- constraints(design.settings)[!grepl(
                            "gc_clamp", names(constraints(design.settings)))]

## ----change_settings2, message = FALSE, warning = FALSE------------------
constraints(design.settings)[["primer_length"]] <- c("min" = 25, "max" = 25)

## ----change_settings3, message = FALSE, warning = FALSE------------------
conOptions(design.settings)[["allowed_mismatches"]] <- 0

## ----write_settings, message = FALSE, warning = FALSE, eval = FALSE------
#  out.file <- tempfile("settings", fileext = ".xml")
#  write_settings(design.settings, out.file)

## ----design_primers, message = FALSE, warning = FALSE, eval = FALSE------
#  optimal.primers <- design_primers(template.df[1:2,], mode.directionality = "fw",
#                                    settings = design.settings)

## ----design_primers_custom, message = FALSE, warning = FALSE, eval = FALSE----
#  optimal.primers.custom <- design_primers(template.df[1:2,], mode.directionality = "fw",
#                                    settings = design.settings, init.algo = "tree",
#                                    opti.algo = "ILP", required.cvg = 0)

## ----store_primers, message = FALSE, warning = FALSE, eval = FALSE-------
#  out.file <- tempfile("my_primers", fileext = ".fasta")
#  write_primers(optimal.primers$opti, out.file)

## ----read_primers, message = FALSE, warning = FALSE----------------------
primer.location <- system.file("extdata", "IMGT_data", "primers", "IGHV", 
                               "Ippolito2012.fasta", package = "openPrimeR")
primer.df <- read_primers(primer.location, fw.id = "_fw", rev.id = "_rev", 
                          merge.ambig = "none", max.degen = 16)

## ----check_constraints, message = FALSE, warning = FALSE-----------------
conOptions(settings)$allowed_other_binding_ratio <- c("max" = 1.0)
constraint.df <- check_constraints(primer.df, template.df, 
                 settings, active.constraints = names(constraints(settings)))

## ----view_primer_coverage, message = FALSE, warning = FALSE--------------
constraint.df$primer_coverage

## ----update_template_cvg, message = FALSE, warning = FALSE---------------
template.df <- update_template_cvg(template.df, constraint.df)

## ----template_coverage, message = FALSE, warning = FALSE-----------------
template.df$primer_coverage[1:5]

## ----cvg_ratio, message = FALSE, warning = FALSE-------------------------
as.numeric(get_cvg_ratio(constraint.df, template.df))

## ----cvg_stats, message = FALSE, warning = FALSE-------------------------
cvg.stats <- get_cvg_stats(constraint.df, template.df, for.viewing = TRUE)

## ----cvg_table, echo=FALSE, results='asis'-------------------------------
knitr::kable(cvg.stats[, !(grepl("_fw", colnames(cvg.stats)) | grepl("_rev", colnames(cvg.stats)))], row.names = FALSE)

## ----template_cvg_plot, fig.show='hold', fig.width=5, fig.height=5-------
plot_template_cvg(constraint.df, template.df)

## ----primer_subsets, message = FALSE, warning = FALSE--------------------
primer.subsets <- subset_primer_set(constraint.df, template.df)

## ----cvg_subsets, fig.show='hold', fig.width=5, fig.height=5-------------
plot_primer_subsets(primer.subsets, template.df)

## ----optimal_subsets, message = FALSE, warning = FALSE-------------------
my.primer.subset <- primer.subsets[[3]]

## ----binding_regions, message = FALSE, warning = FALSE, fig.show='hold', fig.width=5, fig.height=5----
plot_primer_binding_regions(constraint.df, template.df)

## ----cvg_primer_view, message = FALSE, warning = FALSE, fig.show='hold', fig.width=8, fig.height=8----
plot_primer(constraint.df[1,], template.df[1:15,], relation = "fw")

## ----plot_constraints, fig.show='hold', fig.width=7, fig.height=7, message = FALSE, warning = FALSE----
plot_constraint_fulfillment(constraint.df, settings)

## ----constraint_eval_gc, message = FALSE, warning = FALSE----------------
constraint.df$gc_clamp_fw[!constraint.df$EVAL_gc_clamp]
constraints(settings)$gc_clamp

## ----plot_constraint_qualitative, fig.show='hold', fig.width=5, fig.height=5, message = FALSE, warning = FALSE----
plot_constraint(constraint.df, settings, "gc_clamp")

## ----primer_filtering, message = FALSE, warning = FALSE------------------
filtered.df <- filter_primers(constraint.df, template.df, settings,
               active.constraints = c("gc_clamp", "gc_ratio"))

## ----eval_report, message = FALSE, warning = FALSE, eval = FALSE---------
#  my.file <- tempfile()
#  create_report(constraint.df, template.df, my.file,
#                settings, sample.name = "My analysis")

## ----writing_comparison_data, message = FALSE, warning = FALSE-----------
primer.xml <- tempfile("my_primers", fileext =".xml")
write.csv(constraint.df, file = primer.xml, row.names = FALSE)
template.xml <- tempfile("my_templates", fileext = ".xml")
write.csv(constraint.df, file = template.xml, row.names = FALSE)

## ----loading_comparison_data, message = FALSE, warning = FALSE-----------
sel.sets <- c("Glas1999", "Rubinstein1998", "Cardona1995", "Persson1991", "Ippolito2012", "Scheid2011")
primer.files <- list.files(path = system.file("extdata", "IMGT_data", "comparison", 
                           "primer_sets", "IGH", package = "openPrimeR"),
                pattern = "*\\.csv", full.names = TRUE)
primer.data <- read_primers(primer.files)
sel.idx <- which(names(primer.data) %in% sel.sets)
primer.data <- primer.data[sel.idx]
template.files <- rep(system.file("extdata", "IMGT_data", "comparison", "templates", 
                              "IGH_templates.csv", package = "openPrimeR"), 
                              length(primer.data))
template.data <- read_templates(template.files)

## ----comparison_plots_overview, fig.show='hold', fig.width=7, fig.height=7, message = FALSE----
plot_constraint_fulfillment(primer.data, settings, plot.p.vals = FALSE)

## ----comparison_plots_details, fig.show='hold', fig.width=7, fig.height=7, message = FALSE----
plot_constraint(primer.data, settings, active.constraints = c("gc_ratio", "melting_temp_range"))

## ----comparison_primer_binding, fig.show='hold', fig.width=7, fig.height=7, message = FALSE----
plot_primer_binding_regions(primer.data, template.data)


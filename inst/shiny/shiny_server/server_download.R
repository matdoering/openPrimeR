########
# Shiny server functionalities for downloads
#########

observeEvent(input$downloadDataSet, {
	# show a modal informing about missing report dependencies 
	if (input$downloadDataSet == "eval_report" || input$downloadDataSet == "comparison_report") {
		# check rmarkdown render dependencies
		deps <- openPrimeR:::check_report_deps()
		if (any(!deps)) {
			# any depdency is missing
			toggleModal(session, "MissingReportDeps")
		}
	}
})

downloadDataSelector <- reactive({
    # prepare data for download in UI
    # get primers/templates:
    primer.df.raw <- rv_primers$all
    if (length(rv_primers$evaluated_primers) != 0) {
        primer.df.raw <- rv_primers$evaluated_primers
    }
    #template.df <- switch(input$primer_type_download,
        #"all" = rv_templates$cvg_all,
        #"filtered" = rv_templates$cvg_filtered,
        #"optimized" = rv_templates$cvg_optimized)
    # Download xml settings here ...
    # templates may not be annotated with cvg here (not differentiating the type of templates here)
    # n.b.: removed myCatch() as expr() can't be evaluated with try if it's reactive ..
    # n.b.: need to have switch such that 'all' is not always loaded, but only if necessary ..
    # n.b.: `Report` is not part of `all` because it's no data object but directly stored to file.
    if (input$download_style == "raw") {
        # raw output
         primers <- NA
         if (input$set_meta_selector == "all") {
            primers <- primer.df.raw
         } else if (input$set_meta_selector == "filtered") {
            primers <- rv_primers$filtered
         } else {
            primers <- rv_primers$optimized
         }
         out <- switch(input$downloadDataSet,
         ############ ALL ############
                "all" = list(
                       "sequences" = update.sample.name(rv_templates$raw_seqs, input$sample_name),
                       "primers" = primers,
                       "settings" = current.settings()
                       ),
        ########## INDIVIDUAL ############
               "eval_report" = "eval_report",
               "comparison_report" = "comparison_report",
               "coverage_spreadsheet" = "coverage_spreadsheet", 
               "sequences" = update.sample.name(rv_templates$raw_seqs, input$sample_name),
               "primers" = primers,
               "primers_subset" = try(primer_subset_out()),
               "settings" = current.settings(),
               "mismatch_table_fw" = try(mismatch.table.fw()),
               "mismatch_table_rev" = try(mismatch.table.rev()),
               "self_dimerization" = try(self.dimer.data()),
               "cross_dimerization" = try(cross.dimer.data()))
    } else { 
        # formatted output: no downloadinfo added
         primers <- NA
         if (input$set_meta_selector == "all") {
            primers <- primer.df.raw
         } else if (input$set_meta_selector == "filtered") {
            primers <- rv_primers$filtered
         } else {
            primers <- rv_primers$optimized
         }
         # earlier: no FASTA output -> used the tabs output, now i need the raw output ..
         #if (input$set_meta_selector == "all") {
         #   primers <- rv_primers$PrimerTab
         #} else if (input$set_meta_selector == "filtered") {
         #   primers <- rv_primers$PrimerTabFiltered
         #} else {
         #   primers <- rv_primers$PrimerTabOptimized
         #}
         out <- switch(input$downloadDataSet,
            ########### ALL ###############
            "all" = list("sequences" = rv_templates$SeqTab,
                    "primers" = primers,
                    "settings" = current.settings()
                    ),
            ########### INDIVIDUAL ##########
            "eval_report" = "eval_report",
            "comparison_report" = "comparison_report",
            "coverage_spreadsheet" = "coverage_spreadsheet", 
            "sequences" = update.sample.name(rv_templates$raw_seqs, input$sample_name),
            "primers" = primers,
            "primers_subset" = try(view.subset.primers(primer_subset_out(), current.seqs(), run.mode(), input$view_cvg_individual)),
            "settings" = current.settings(),
            "mismatch_table_fw" = try(view.mismatch.table(mismatch.table.fw())),
            "mismatch_table_rev" = try(view.mismatch.table(mismatch.table.rev())), 
            "self_dimerization" = try(self.dimer.data()),
            "cross_dimerization" = try(cross.dimer.data()))
        }
    return(out)
})

output$downloadData <- downloadHandler(
    # create download data depending on data selected in downloadDataSelector()
    filename = function() { 
        # return file where we want to save
        data <- downloadDataSelector()
        out <-  paste("openPrimeR_", input$sample_name, "_", input$downloadDataSet, sep = "")
        if (length(data) == 0) {
            out <- paste(out, ".error", sep = "")
        } else if (is(data, "DesignSettings")) {
            out <- paste(out, ".xml", sep = "")
        } else if (unlist(data)[1] %in% c("eval_report", "comparison_report")) {
            # pdf report
            out <- paste(out, ".pdf", sep = "")
        } else if (unlist(data)[1] %in% c("coverage_spreadsheet")) {
            # xls output
            out <- paste(out, ".xlsx", sep = "")
        } else if (class(data) == "list") {
            out <- paste(out, ".zip", sep = "")
        } else if ((is(data, "Primers") || is(data, "Templates")) && input$download_style != "raw") {
            # FASTA output
            out <- paste(out, ".fasta", sep = "")
        } else if (length(nrow(data)) != 0) { # csv output
            out <- paste(out, ".csv", sep = "")
        } else {
            # unknown data type
            out <- paste(out, ".txt", sep = "")
        }
        return(out)
    },
    content = function(file) { # write output to file
        data <- downloadDataSelector()
        #print("input data")
        #print(data)
        if (length(data) == 0) {
            write("", file)
            return()
        }
        # get primers/templates:
        primer.df <- switch(input$set_meta_selector,
            "all" = rv_primers$evaluated_primers,
            "filtered" = current.filtered.primers(),
            "optimized" = optimal.primers())
        template.df <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
        if (is(data, "DesignSettings")) {
            openPrimeR:::write_settings(data, file)
        } else if (is(data, "Primers") && input$download_style != "raw") {
            # write FASTA output
            openPrimeR:::write_primers(data, file)
        } else if (is(data, "Templates") && input$download_style != "raw") {
            # write FASTA output
            openPrimeR:::write_templates(data, file)
        } else if (is(data, "Primers") && input$download_style == "raw") {
            # write CSV output
            openPrimeR:::write_primers(data, file, "CSV")
        } else if (is(data, "Templates") && input$download_style == "raw") {
            # write CSV output
            openPrimeR:::write_templates(data, file, "CSV")
        } else if (unlist(data)[1] == "eval_report") {
            if (rmarkdown::pandoc_available()) {
                res <- try(openPrimeR::create_report(primer.df, template.df, 
                            file, current.settings(), 
                            sample.name = input$sample_name, 
                            used.settings = rv_primers$optimal_data$used_constraints))
            } 
        } else if (unlist(data)[1] == "comparison_report") {
            if (rmarkdown::pandoc_available()) {
                res <- try(openPrimeR::create_report(
                            rv_comparison.data$comp_primers, 
                            rv_comparison.data$comp_templates, 
                            file, current.settings(), 
                            sample.name = input$sample_name))
            } else {
                msg <- paste0("Pandoc for rmarkdown is not available on your system.", 
                        "Please install it first to generate a report.")
                warning(msg)
            }
        } else if (unlist(data)[1] == "coverage_spreadsheet") {
            res <- try(openPrimeR::create_coverage_xls(
                            primer.df, template.df, file,
                            current.settings()))
        } else if (class(data) == "list") {
            # "all" output -> zip
            tmpdir <- tempdir() 
            fs <- rep(NA, length(data))
            for (i in seq_along(data)) {
                if (length(data[[i]]) != 0) { # non-null data
                    if (is(data[[i]], "DesignSettings")) {
                        # settings object
                        fs[i] <- file.path(tmpdir, paste(names(data)[i], ".xml", sep = ""))
                        openPrimeR:::write_settings(data[[i]], fs[i])
                    } else if (length(nrow(data[[i]])) == 0 && length(data[[i]]) != 0) {
                        # other type of data
                        fs[i] <- file.path(tmpdir, paste(names(data)[i], ".txt", sep = ""))
                        write(data[[i]], fs[i])
                    } else if (nrow(data[[i]]) != 0) { 
                        # data frame -> write csv
                        fs[i] <- file.path(tmpdir, paste(names(data)[i], ".csv", sep = ""))
                        write.csv(data[[i]], fs[i], row.names = FALSE)
                    } 
                }
            }
            zip(zipfile=file, files=fs, flags = "-q -j")
            unlink(fs)
        } else if (length(nrow(data)) == 0) {
            # text output
            write(data, file)
        } else if (nrow(data) != 0) {
            # data frame
            write.csv(data, file, row.names = FALSE)
        }  
        # finally, check whether any file was created.
        if (!file.exists(file)) {
            # something has gone wrong (no data/error)
            write("", file)
        }
})

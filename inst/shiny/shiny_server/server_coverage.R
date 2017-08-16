#########
# Shiny server functionalities relating to primer coverage
##########

primerViewObserverGroup <- observeEvent(current.seqs(), {
    # update coverage group selector when templates are available
    if (length(current.seqs()) != 0) {
        groups <- unique(current.seqs()$Group)
        updateSelectInput(session, "selected_group_coverage", choices = c("all", groups))
    }
})

primerViewObserver <- observe({
    # update input$selected_primer for coverage tab
    cur.table <- NULL
    if (input$set_meta_selector == "all") {
        if (length(rv_primers$evaluated_primers) == 0) {
            updateSelectInput(session, "selected_primer", choices = "")
            return()
        }
        cur.table <- rv_primers$evaluated_primers
    } else if (input$set_meta_selector == "filtered") {
        cur.table <- current.filtered.primers()
        if (length(cur.table) == 0) {
            updateSelectInput(session, "selected_primer", choices = "")
            return() 
        }
    } else if (input$set_meta_selector == "optimized") {
        cur.table <- optimal.primers()
        if (length(cur.table) == 0) {
            updateSelectInput(session, "selected_primer", choices = "")
            validate(need(FALSE, "Please compute the optimized data set first."))
            return()
        }
    } else {
        updateSelectInput(session, "selected_primer", choices = "")
        validate(need(FALSE, "No valid primer set selected."))
        return()
    }
    template.df <- current.seqs()
    if (length(cur.table) != 0 && nrow(cur.table) != 0 && "primer_coverage" %in% colnames(cur.table)) {
        cov.seqs <- sapply(strsplit(cur.table$Covered_Seqs, split = ","), function(x) match(as.numeric(x), template.df$Identifier))
        groups <- sapply(seq_along(cov.seqs), function(x) unique(template.df$Group[cov.seqs[[x]]]))
        if (length(input$selected_group_coverage) == 0 || any(c("", "all") %in% input$selected_group_coverage)) {
            sel <- seq_len(nrow(cur.table))
        } else {
            sel <- which(sapply(groups, function(x) input$selected_group_coverage %in% x))
        }
        cur.table <- cur.table[sel,]
        updateSelectInput(session, "selected_primer", choices =  cur.table$ID)
    }
})

coverage.statistics <- reactive({
    # overview of coverage information
    data <- switch(input$set_meta_selector,
            "all" = rv_primers$evaluated_primers,
            "filtered" = current.filtered.primers(),
            "optimized" = optimal.primers())
    template.data <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    validate(need(data, "No primer coverage available."))
    validate(need(template.data, "No templates available."))
    # stats for expected coverage:
    cvg.stats.exp <- openPrimeR:::get_cvg_stats(data, template.data)
    if (length(cvg.stats.exp) != 0 && nrow(cvg.stats.exp) != 0) {
        cvg.stats.exp <- cbind(cvg.stats.exp, Coverage_Definition = "Expected_Coverage")
    }
    # stats for text identity coverage:
    cvg.stats.txt <- openPrimeR:::get_cvg_stats(data, template.data,
                        allowed.mismatches = 0, cvg.definition = "basic")
    if (length(cvg.stats.txt) != 0 && nrow(cvg.stats.txt) != 0) {
        cvg.stats.txt <- cbind(cvg.stats.txt, Coverage_Definition = "Identity_Coverage")
    }
    cvg.stats <-  rbind(cvg.stats.exp, cvg.stats.txt)
    return(cvg.stats)
})
coverage.statistics.mismatches <- reactive({
    # overview of coverage information for different mismatch settings
    data <- switch(input$set_meta_selector,
            "all" = rv_primers$evaluated_primers,
            "filtered" = current.filtered.primers(),
            "optimized" = optimal.primers())
    template.data <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    validate(need(data, "No primer coverage available."))
    validate(need(template.data, "No templates available."))
    # stats for expected coverage:
    mm.range <- seq(0, input$allowed_mismatches)
    mm.stats <- vector("list", length(mm.range))
    for (i in seq_along(mm.range)) {
        cvg.stats <- openPrimeR:::get_cvg_stats(data, template.data,
                        allowed.mismatches = mm.range[i])
        cvg.stats <- cbind(cvg.stats, "MaxMismatches" = mm.range[i])
        mm.stats[[i]] <- cvg.stats
    }
    cvg.string <- sapply(mm.stats, function(x) paste0(round(x[x$Group == "Total", "Coverage_Ratio"] * 100, 1), "%"))
    mm.stats <- do.call(rbind, mm.stats)
    # update the selector for max mismatches:
    labels <- paste0(mm.range, " mismatches (", cvg.string, " coverage)")
    opts <- mm.range
    names(opts) <- labels
    updateSelectInput(session, "allowed_mm_cvg_stats", choices = opts)
    return(mm.stats)
})

#basic_string_cvg <- reactive({
    ## basic coverage according to string complementarity with 0 mismatches
    #data <- switch(input$set_meta_selector,
            #"all" = rv_primers$evaluated_primers,
            #"filtered" = current.filtered.primers(),
            #"optimized" = optimal.primers())
    #template.data <- switch(input$set_meta_selector,
            #"all" = rv_templates$cvg_all,
            #"filtered" = rv_templates$cvg_filtered,
            #"optimized" = rv_templates$cvg_optimized)
    #validate(need(data, "No primer coverage available."))
    #validate(need(template.data, "No templates available."))
    #cvg <- openPrimeR::get_cvg_ratio(data, template.data, allowed.mismatches = 0, cvg.definition = "basic")
    #return(cvg)
#})

output$cvg_stats <- DT::renderDataTable({
    # table with coverage stats
    data <- switch(input$set_meta_selector,
            "all" = rv_primers$evaluated_primers,
            "filtered" = current.filtered.primers(),
            "optimized" = optimal.primers())
    template.data <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    validate(need(data, "No primer coverage available."))
    validate(need(template.data, "No templates available."))
    # select the coverage definition for template coverage table:
    selected.cvg.def <- input$selected_cvg_def_stats
    stats <- coverage.statistics()
    if (length(stats) != 0) {
        stats <- stats[stats$Coverage_Definition == selected.cvg.def,]
        stats <- stats[, c("Group", "Coverage", "Coverage_fw", "Coverage_rev")]
    }
    return(DT::datatable(stats,
        rownames = FALSE, options = list(dom = "pt"), # dom = pt -> show pages and table
        caption = paste("The number of covered sequences per group of templates.",
            "If the coverage definition is set to 'Expected Coverage', ",
            "coverage is computed using the extended coverage criteria",
            "and (possibly) allowing for multiple mismatches.",
            "If, however, the coverage definition is set to",
            "'Identity Coverage', only coverage events where the primers",
            "are perfectly complementary to the templates are considered.")))
})
output$cvg_stats_mismatches <- DT::renderDataTable({
    # table with template coverage stats
    data <- switch(input$set_meta_selector,
            "all" = rv_primers$evaluated_primers,
            "filtered" = current.filtered.primers(),
            "optimized" = optimal.primers())
    template.data <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    validate(need(data, "No primer coverage available."))
    validate(need(template.data, "No templates available."))
    # select the coverage definition for template coverage table:
    allowed.mm <- input$allowed_mm_cvg_stats
    stats <- coverage.statistics.mismatches() # TODO: this is new
    if (allowed.mm == "") {
        # nothing to select
        return(NULL)
    } else {
        # convert from string to numeric
        allowed.mm <- as.numeric(allowed.mm)
    }
    if (length(stats) != 0) {
        stats <- stats[stats$MaxMismatches == allowed.mm,]
        stats <- stats[, c("Group", "Coverage", "Coverage_fw", "Coverage_rev")]
    }
    return(DT::datatable(stats,
        rownames = FALSE, options = list(dom = "pt"), # dom = pt -> show pages and table
        caption = paste("The number of sequences that are expected to be",
                    "covered by primers when allowing for a certain number",
                    "of mismatches.")))
})


output$CoverageTotal <- renderUI({
    # output text giving an overview of coverage/Tm Delta

    if (length(coverage.statistics()) == 0) {
        return(NULL)
    }
    stats <- coverage.statistics()
    selected.group <- input$selected_group_coverage
    if (length(stats) != 0) {
        cvg.stats.exp <- stats[stats$Coverage_Definition == "Expected_Coverage",]
        cvg.stats.ident <- stats[stats$Coverage_Definition == "Identity_Coverage",]
        cvg.text.ident <- paste("<li>", 
                        openPrimeR:::create.cvg.text(cvg.stats.ident, selected.group, "Identity Coverage"), 
                        "</li>", sep = "")
        cvg.text.exp <- paste("<li>", 
                        openPrimeR:::create.cvg.text(cvg.stats.exp, selected.group, "Expected Coverage"), 
                        "</li>", sep = "")

    } else {
        cvg.text.ident <- ""
        cvg.text.exp <- ""
    }

    text <- HTML(paste("<h3><ul>", 
                 cvg.text.exp,
                 cvg.text.ident,
                 "</ul></h3>", sep = ""))
    return(text)
})
output$ConstraintsTotal <- renderUI({
    Tm.info <- cur.Tm.info()
    if (length(Tm.info) != 0) {
        Tm.text <- paste("<li>", create.Tm.text(Tm.info$Tm_range, Tm.info$Tm_diff), "</li>", sep = "")
    } else {
        Tm.text <- ""
    }
    annealing.text <- ""
    if (length(annealing.temperature()) != 0) {
        opt <- ifelse(input$automatic_annealing_temp == "active", "recommended", "user-defined")
        annealing.text <- paste0("<li>The ", opt, " annealing temperature is ", annealing.temperature(), "&#8451;.</li>")
    }
    text <- HTML(paste("<h3><ul>", 
                 Tm.text,
                 annealing.text,
                 "</ul></h3>", sep = ""))
    return(text)
})

primer_plot_height <- reactive({
    # height of the primer view plot
    primer.df <- switch(input$set_meta_selector, 
                    "all" = rv_primers$evaluated_primers,
                    "filtered" = current.filtered.primers(),
                    "optimized" = optimal.primers())
    template.df <- current.seqs()
    if (length(primer.df) == 0 || length(template.df) == 0) {
        return(0)
    }
    id <- NULL
    if (length(input$selected_primer) == 0 || input$selected_primer == "") { # no ID given -> plot all
        id <- primer.df$Identifier
    } else {
        id <- primer.df$Identifier[match(input$selected_primer, primer.df$ID)] # match from primer identifier to ID
    }
    m <- match(id, primer.df$Identifier)
    primer.df <- primer.df[m, ]
    show.group <- input$selected_group_coverage
    if (!is.null(show.group) && !"all" %in% show.group) {
        lex.id <- which(template.df$Group %in% show.group)
        template.df <- template.df[lex.id,]
    }
    if ("primer_coverage" %in% colnames(primer.df)) {
        n1 <- sum(primer.df$primer_coverage) # nbr of coverage events to show
    } else {
        n1 <- 0
    }
    n2 <- nrow(template.df) # nbr of templates to show
    n <- n1 + n2
    height <- openPrimeR:::get.plot.height(n, px.per.n = 30)
    return(height)
})

primer_plot_width <- reactive({
    # width of the primer view plot
    myData <- current.seqs()
    if (length(myData) == 0) {
        return(1200)
    }
    n <- max(nchar(myData$Sequence))
    width <- openPrimeR:::get.plot.height(n, px.per.n = 2)
    return(width)
})


output$primer_plot <- renderPlot({
    # data for the primer view plot
    primer.df <- switch(input$set_meta_selector, 
                    "all" = rv_primers$evaluated_primers,
                    "filtered" = current.filtered.primers(),
                    "optimized" = optimal.primers())
    validate(need(primer.df, "The selected primer data set cannot be plotted, because no primer coverage is available in the selected set."))
    validate(need(nrow(primer.df) !=0, "The selected primer data set cannot be plotted, because no primers are available."))
    validate(need(current.seqs(), "No template sequences are available for plotting yet."))
    show.group <- input$selected_group_coverage
    template.df <- current.seqs()
    id <- NULL
    if (length(input$selected_primer) == 0 || input$selected_primer == "") { # no ID given -> plot all
        id <- primer.df$Identifier
    } else {
        id <- primer.df$Identifier[match(input$selected_primer, primer.df$ID)] # match from primer identifier to ID
    }
    if (!is.null(show.group) && !"all" %in% show.group) {
        lex.id <- which(template.df$Group %in% show.group)
        template.df <- template.df[lex.id,]
    }
    validate(need(template.df, "No template sequences avaialable."))
    # check whether selected primers have any cvg?
    m <- match(id, primer.df$Identifier)
    validate(need(any(primer.df[m, "primer_coverage"] != 0), "Selected primers do not cover any template sequences."))
    selected.group <- "all" # was 'input$selected_group' before
    relation <- input$primer_plot_rel
    region.names <- NULL
    # change region names if we have immuno data
    if (input$template_scenario == "supplied" && input$selected_supplied_templates == "immunological") {
        openPrimeR:::plot_primer(primer.df, template.df, id, relation, region.names = c("Leader", "Variable region"))
    } else {
        openPrimeR:::plot_primer(primer.df, template.df, id, relation)
    }
},  width = primer_plot_width, height=primer_plot_height, units="px")

cvg.group.plot.dim <- reactive({
    # dimension of group cvg plot
    primer.df <- switch(input$set_meta_selector,
        "all" = rv_primers$evaluated_primers,
        "filtered" = current.filtered.primers(),
        "optimized" = optimal.primers())
    template.df <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    if (length(template.df) == 0 || length(primer.df) == 0)  {
        #print("templates unknown for group plot dimension ...")
        return(list("width" = 800, "height" = 800))
    }
    lex.sel <- input$selected_group_coverage
    if (!is.null(lex.sel) && !"all" %in% lex.sel) { # select subset
        idx <- which(template.df$Group %in% lex.sel)
        template.df <- template.df[idx,]
    }
    nbr.groups <- length(unique(template.df$Group))
    width <- openPrimeR:::get.plot.height(nbr.groups, 25, 500) # nbr of gruops
    max.nbr.templates.per.group <- max(table(template.df$Group))
    height <- openPrimeR:::get.plot.height(max.nbr.templates.per.group, 1, 500)
    out <- list("width" = width, "height" = height)
    return(out)
})

cvg.group.plot.width <- reactive({
    # width for group coverage plot
    if (length(cvg.group.plot.dim()) == 0) {
        return(NULL)
    }
    return(cvg.group.plot.dim()$width)
})
cvg.group.plot.height <- reactive({
    # height for group coverage plot
    if (length(cvg.group.plot.dim()) == 0) {
        return(NULL)
    }
    return(cvg.group.plot.dim()$height)
})

output$Coverage_Group <- renderPlot({
    # group coverage plot
    lex.sel <- input$selected_group_coverage
    primer.data <- switch(input$set_meta_selector,
            "all" = rv_primers$evaluated_primers,
            "filtered" = current.filtered.primers(),
            "optimized" = optimal.primers())
    template.data <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    validate(need(primer.data, "No primer coverage available."))
    validate(need(template.data, "No templates available."))
    openPrimeR:::plot_template_cvg(primer.data, template.data, groups = lex.sel)
}, width = cvg.group.plot.width, height = cvg.group.plot.height)

output$coverage_primer_per_group_ui <- renderUI({
    # ui output of group coverage; prevents overplotting for multiple plot elements on one page
    plotOutput("Coverage_Group",width = paste0(cvg.group.plot.width(), "px"), height = paste0(cvg.group.plot.height(), "px"))
})
cvg.stats.primer.mismatches <- reactive({
    data <- switch(input$set_meta_selector,
            "all" = rv_primers$evaluated_primers,
            "filtered" = current.filtered.primers(),
            "optimized" = optimal.primers())
    template.data <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    validate(need(data, "No primer coverage available."))
    validate(need(template.data, "No templates available."))
    cvg.stats <- openPrimeR:::get_cvg_stats_primer(data, template.data)
    return(cvg.stats)
})
output$primer_cvg_stats_mismatch <- DT::renderDataTable({
    # table with coverage stats wrt mismatches
    primer.df <- switch(input$set_meta_selector,
            "all" = rv_primers$evaluated_primers,
            "filtered" = current.filtered.primers(),
            "optimized" = optimal.primers())
    template.df <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    validate(need(primer.df, "No primer coverage available."))
    validate(need(template.df, "No templates available."))
    return(DT::datatable(cvg.stats.primer.mismatches(),
           caption = paste("The number of covered template sequence for every primer.",
                            "Columns with numeric identifiers give the number of coverage events",
                            "corresponding to primers binding with a certain number of mismatches.",
                            "The group coverage indicates the percentage of covered templates",
                            "from each group."),
           rownames = FALSE, options = list(dom = "pt"))) # dom = pt -> show pages and table
})

output$template_coverage_mismatch_ui <- renderUI({
    # important: set size of plot here to prevent overlap in the UI
    plotOutput("template_coverage_mismatch", width = 1200, height = cvg.template.mismatch.height()) # TODO: change height!!
})
cvg.template.mismatch.height <- reactive({
    # height of template coverage plot, stratified by mismatches
    primer.df <- switch(input$set_meta_selector,
        "all" = rv_primers$evaluated_primers,
        "filtered" = current.filtered.primers(),
        "optimized" = optimal.primers())
    template.df <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    if (length(primer.df) == 0 || length(template.df) == 0 || !"primer_coverage" %in% colnames(primer.df))  {
        return(1200)
    }
    max.cvg <- openPrimeR:::get_cvg_ratio(primer.df, template.df)
    max.mm <- max(as.numeric(unlist(strsplit(c(primer.df$Nbr_of_mismatches_fw, primer.df$Nbr_of_mismatches_rev), split = ","))))
    for (i in seq(0, max.mm - 1)) {
        cur.cvg <- openPrimeR:::get_cvg_ratio(primer.df, template.df, allowed.mismatches = i)
        if (cur.cvg == max.cvg) {
            max.mm <- i
            break
        }
    }
    height <- openPrimeR:::get.plot.height(ceiling(max.mm + 1) / 2, 400, 600)
    return(height)
})

output$template_coverage_mismatch <- renderPlot({
    # plot of template coverage with consideration of mismatch setting
     primer.df <- switch(input$set_meta_selector,
        "all" = rv_primers$evaluated_primers,
        "filtered" = current.filtered.primers(),
        "optimized" = optimal.primers())
    validate(need(primer.df$primer_coverage, "Primer coverage has not been computed yet."))
    validate(need(any(primer.df$primer_coverage != 0), "The primers do not cover any of the templates."))
    template.df <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    validate(need(template.df, "The template data set is not available."))
    lex.sel <- input$selected_group_coverage
    if (!is.null(lex.sel) && !"all" %in% lex.sel) {
        # select a template subset
        idx <- which(template.df$Group %in% lex.sel)
        template.df <- template.df[idx,]
    }
    openPrimeR:::plot_template_cvg(primer.df, template.df, per.mismatch = TRUE)
})


cvg.primer.group.width <- reactive({
    # width of individual primer plot
     primer.df <- switch(input$set_meta_selector,
        "all" = rv_primers$evaluated_primers,
        "filtered" = current.filtered.primers(),
        "optimized" = optimal.primers())
    if (length(primer.df) == 0)  {
        print("primers unknown for group plot dimension ...")
        return(800)
    }
    width <- openPrimeR:::get.plot.height(nrow(primer.df), 20, 500)
    return(width)
})
cvg.primer.group.height <- reactive({
    # height of individual primer plot
    return(600)
})
output$Coverage_Primer_mismatches <- renderPlot({
    # individual primer coverage plot: mismatches
    primer.df <- switch(input$set_meta_selector,
        "all" = rv_primers$evaluated_primers,
        "filtered" = current.filtered.primers(),
        "optimized" = optimal.primers())
    validate(need(primer.df$primer_coverage, "Primer coverage has not been computed yet."))
    validate(need(any(primer.df$primer_coverage != 0), "The primers do not cover any of the templates."))
    template.df <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    validate(need(template.df, "The template data set is not available."))
    lex.sel <- input$selected_group_coverage
    if (!is.null(lex.sel) && !"all" %in% lex.sel) { # select subset
        idx <- which(template.df$Group %in% lex.sel)
        template.df <- template.df[idx,]
    }
    openPrimeR:::plot_primer_cvg(primer.df, template.df, per.mismatch = TRUE)
}, width = 1200, height = 1200)

output$Coverage_Primer <- renderPlot({
    # individual primer coverage plot
    primer.df <- switch(input$set_meta_selector,
        "all" = rv_primers$evaluated_primers,
        "filtered" = current.filtered.primers(),
        "optimized" = optimal.primers())
    validate(need(primer.df$primer_coverage, "Primer coverage has not been computed yet."))
    validate(need(any(primer.df$primer_coverage != 0), "The primers do not cover any of the templates."))
    template.df <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    validate(need(template.df, "The template data set is not available."))
    openPrimeR:::plot_primer_cvg(primer.df, template.df, groups = input$selected_group_coverage)
})#, width = cvg.primer.group.width, height = cvg.primer.group.height)

output$coverage_primer_groups_ui <- renderUI({
    # ui output of individual primer coverage plot
    plotOutput("Coverage_Primer")#, width = paste0(cvg.primer.group.width(), "px"), height = paste0(cvg.primer.group.height(), "px"))
})


output$primer_binding_regions <- renderPlot({
    # plot of primer binding regions
    data <- switch(input$set_meta_selector,
        "all" = rv_primers$evaluated_primers,
        "filtered" = current.filtered.primers(),
        "optimized" = optimal.primers())
    validate(need(data, "Please compute the primer coverage first."))
    template.df <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    validate(need(template.df, "Please make sure the templates are present."))
    lex.sel <- input$selected_group_coverage
    if (length(template.df) == 0 || length(data) == 0) {
        return(NULL)
    }
    # why select a subset beforehand?
    relation <- input$primer_location_plot_rel
    # change annotation of regions for immuno data 
    if (input$template_scenario == "supplied" && input$selected_supplied_templates == "immunological") {
        openPrimeR:::plot_primer_binding_regions(data, template.df, relation = relation, group = lex.sel, 
                                                region.names = c("Leader", "Variable region"))
    } else {
        openPrimeR:::plot_primer_binding_regions(data, template.df, relation = relation, group = lex.sel)
    }
})

mismatch.table.fw <- reactive({
    # table with mismatch data for fw primers
    if (input$set_meta_selector == "all") {
            primer.df <- rv_primers$evaluated_primers
            template.df <- rv_templates$cvg_all
        } else if (input$set_meta_selector == "filtered") {
            primer.df <- current.filtered.primers()
            template.df <- rv_templates$cvg_filtered
        } else if (input$set_meta_selector == "optimized") {
            primer.df <- optimal.primers()
            template.df <- rv_templates$cvg_optimized
        } else {
            return(NULL)
        }
        validate(need(primer.df$primer_coverage, "Mismatch data (fw) is not available since primer coverage has not been computed yet."))
        validate(need(template.df, "Mismatch data (fw) is not available since no templates are available."))
        table <- openPrimeR:::compute.mismatch.table(primer.df, template.df, "fw")
        validate(need(table, "No primers (fw) binding with mismatches found.")) 
        return(table)
})
mismatch.table.rev <- reactive({
    # table with mismatch data for rev primers
    if (input$set_meta_selector == "all") {
            primer.df <- rv_primers$evaluated_primers
            template.df <- rv_templates$cvg_all
        } else if (input$set_meta_selector == "filtered") {
            primer.df <- current.filtered.primers()
            template.df <- rv_templates$cvg_filtered
        } else if (input$set_meta_selector == "optimized") {
            primer.df <- optimal.primers()
            template.df <- rv_templates$cvg_optimized
        } else {
            return(NULL)
        }
        validate(need(primer.df, "Mismatch data (rev) is not available since primer coverage has not been computed yet."))
        validate(need(template.df, "Mismatch data (rev) is not available since no templates are available."))
        table <- openPrimeR:::compute.mismatch.table(primer.df, template.df, "rev")
        validate(need(table, "No primers (rev) binding with mismatches found."))
        return(table)
})
display.mismatch.table <- reactive({
    # the currently selected mismatch table
        if (input$set_meta_selector == "all") {
            primer.df <- rv_primers$evaluated_primers
            template.df <- rv_templates$cvg_all
        } else if (input$set_meta_selector == "filtered") {
            primer.df <- current.filtered.primers()
            template.df <- rv_templates$cvg_filtered
        } else if (input$set_meta_selector == "optimized") {
            primer.df <- optimal.primers()
            template.df <- rv_templates$cvg_optimized
        } else {
            return(NULL)
        }
        validate(need(primer.df, "No primers available."))
        validate(need(template.df, "No templates available."))
        withProgress(message = 'Rendering mismatch table ...', value = 0, {
            if (input$selected_primer_set_mismatches_direction == "fw") {
                table <- mismatch.table.fw()
            } else {
                table <- mismatch.table.rev()
            }
        })
        return(view.mismatch.table(table))
    })
output$mismatch_table <- DT::renderDataTable(
    # shows the currently selected mismatch table 
    DT::datatable(
        display.mismatch.table(),
        , escape=FALSE, 
        caption = "Overview of all primers binding with mismatches. The input template sequences (top rows) as well as the amplicons resulting from mismatch primers (bottom rows) are displayed in their nucleotide and amino-acid sequence."
    , options = list(), rownames=FALSE
))

cur.Tm.info <- reactive({
    # current maximal difference in primer Tm's
    data <- switch(input$set_meta_selector,
            "all" = rv_primers$evaluated_primers,
            "filtered" = current.filtered.primers(),
            "optimized" = optimal.primers())
    if (length(data) == 0 || !"melting_temp" %in% colnames(data)) {
        return(NULL)
    }
    Tm.range <- c(min(data$melting_temp), max(data$melting_temp))
    max.diff <- Tm.range[2] - Tm.range[1]
    Tm.info <- list("Tm_range" = Tm.range, "Tm_diff" = max.diff)
    return(Tm.info)
})

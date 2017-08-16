######### Shiny server functionalities for filtering

current.filtered.primers <- reactive({
    # current filtered primers: either after evaluation / direct filtering /
    # optimization
    quick.filtered <- rv_primers$filtered$data
    opti.filtered <- rv_primers$filtered_opti$data  # the resulting primers from filtering during the optimization procedure
    val <- rv_primers$last_filter_type  # use 'last_filter_type' to decide which are the most recent filtered primers
    ret <- NULL
    if (length(val) != 0) {
        if (val == "filtered" && length(quick.filtered) != 0) {
            ret <- quick.filtered
        } else if (val == "opti" && length(opti.filtered) != 0) {
            ret <- opti.filtered
        }
    }
    return(ret)
})

filtered.set.up.to.date <- reactive({
    # determine whether primer set was filtered using the current constraints or not
    if (length(rv_values$last_filtering_constraints) != 0) {
        if (!openPrimeR:::compare.constraints(rv_values$last_filtering_constraints, 
            constraints()[["active_settings"]])) {
            return(FALSE)
        } else {
            return(TRUE)
        }
    } else {
        return(FALSE)
    }
})

output$filtering.constraints.status <- renderText({
    # output whether the primer set was filtered using the current constraints or not
    filtered.msg <- "Primer set has been filtered using the selected constraints."
    not.filtered.msg <- "Primer set has NOT been filtered using the selected constraints yet."
    out <- NULL
    if (filtered.set.up.to.date()) {
        out <- filtered.msg
    } else {
        out <- not.filtered.msg
    }
    return(out)
})

filtered.primers.quick <- observeEvent(input$quickFilterButton, {
    # start cascaded filter when quickfilterButton is pressed
    if (input$quickFilterButton == 0) {
        return(NULL)
    }
    if (length(primer.data()) == 0 || length(current.seqs()) == 0) {
        session$sendCustomMessage(type = "jsCode", list(value = "$('#NotifyNoDataAvailable').modal('show')"))
        return(NULL)
    }
    primer.df <- primer.data()
    if (length(rv_primers$evaluated_primers) != 0) {
        # if primers were evaluated -> use the evaluated data set to save some
        # computations
        primer.df <- rv_primers$evaluated_primers
    }
    ######### progress
    progress <- shiny::Progress$new()
    progress$set(message = "Filtering", value = 0)
    on.exit(progress$close())
    updateProgress <- function(value, detail, option) {
        if (is.null(value)) {
            value <- progress$getValue()
            value <- value + (progress$getMax() - value)/5
        }
        if (option == "inc") {
            progress$inc(amount = value, detail = detail)
        } else {
            progress$set(value = value, detail = detail)
        }
    }
    # only compute constraints in filtering that weren't determined already
    active.constraints <- names(openPrimeR:::constraints(current.settings()))
    comp.constraints <- setdiff(active.constraints, isolate(rv_primers$available_constraints))
    filtered.data <- openPrimeR:::cascaded.filter.quick(primer.df, 
        current.seqs(), current.settings(), 
        to.compute.constraints = comp.constraints,
        mode.directionality = run.mode(),
        active.constraints = active.constraints,
        no.structures = FALSE,
        updateProgress = updateProgress)
    rv_templates$cvg_filtered <- openPrimeR:::update_template_cvg(current.seqs(), 
                                                                  filtered.data$data, run.mode())
    rv_primers$last_filter_type <- "filtered"  # indicate the most recent table: filtered data
    # save applied filtering constraints
    rv_values$last_filtering_constraints <- openPrimeR:::filters(current.settings())
    rv_primers$filtered <- filtered.data
    # update tabsets w/ filtering results
    switch.view.selection("filtered", input$main, session)
})

output$filtering_stats_cvg <- renderPlot({
    # plot cvg statistics on the filtering procedure
    if (input$set_meta_selector == "filtered") {
        stats <- rv_primers$filtered$stats
        stats.relax <- NULL
    } else if (input$set_meta_selector == "optimized") {
        validate(need(rv_primers$filtered_opti, "Filtered optimization primers are not available."))
        stats <- rv_primers$filtered_opti$stats
        stats.relax <- rv_primers$filtered_opti$stats_relax
    } else {
        return(NULL)
    }
    validate(need(stats, "No stats available, because no filters were active."))
    validate(need(any(!is.na(stats$Current_Coverage)), "No stats available, because not coverage was computed."))
    openPrimeR:::plot.filtering.stats.cvg(stats, stats.relax)
})

output$filtering_stats <- renderPlot({
    # plot nbr of excluded primers during the filtering procedure
    if (input$set_meta_selector == "filtered") {
        stats <- rv_primers$filtered$stats
        stats.relax <- NULL
    } else if (input$set_meta_selector == "optimized") {
        validate(need(rv_primers$filtered_opti, "Filtered optimization primers are not available."))
        stats <- rv_primers$filtered_opti$stats
        stats.relax <- rv_primers$filtered_opti$stats_relax
    } else {
        return(NULL)
    }
    validate(need(stats, "No stats available, because no filters were active."))
    openPrimeR:::plot.filtering.stats(stats, stats.relax)
})

output$filtering_runtime <- renderPlot({
    # show the runtime of individual filtering steps
    if (input$set_meta_selector == "filtered") {
        stats <- rv_primers$filtered$stats
    } else if (input$set_meta_selector == "optimized") {
        validate(need(rv_primers$filtered_opti, "Filtered optimization primers are not available"))
        stats <- rv_primers$filtered_opti$stats
    } else {
        return(NULL)
    }
    validate(need(stats, "No stats available, because no filters were active."))
    openPrimeR:::plot.filtering.runtime(stats)
})

output$exclusion_stats <- renderPlot({
    # show statistics on excluded primers during filtering
    if (input$set_meta_selector == "filtered") {
        excluded <- rv_primers$filtered$excluded
        stats <- rv_primers$filtered$stats
    } else if (input$set_meta_selector == "optimized") {
        validate(need(rv_primers$filtered_opti, "Filtered optimization primers are not available."))
        excluded <- rv_primers$filtered_opti$excluded
        stats <- rv_primers$filtered_opti$stats
    } else {
        return(NULL)
    }
    validate(need(stats, "No stats available, because no filters were active."))
    validate(need(current.seqs(), "No template sequences available."))
    if (length(excluded) == 0 || length(stats) == 0 || !"primer_coverage" %in% colnames(excluded)) {
        # can't create plot without cvg
        return(NULL)
    }
    openPrimeR:::plot.excluded.hist(excluded, stats, current.seqs())
})

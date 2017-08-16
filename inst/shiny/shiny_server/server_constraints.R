###### Server functionalities for comparing primers
constraint.plot.height <- reactive({
    nbr.rows <- length(input$selected_other_result) / 2 # 2 columns
    width <- openPrimeR:::get.plot.height(nbr.rows, 200, 200)
    return(width)
})
output$constraint_plot_histogram <- renderPlot({
    # plot melting temperatures
    data <- switch(input$set_meta_selector, all = rv_primers$evaluated_primers, 
        filtered = current.filtered.primers(), optimized = optimal.primers())
    validate(need(data, "No data available for plotting."))
    p <- openPrimeR:::plot_constraint(data, current.settings(), input$selected_other_result)
    return(p)
}, height = constraint.plot.height)

output$constraint_plots_cvg_constraints <- renderPlot({
    # plot cvg constraints
    data <- switch(input$set_meta_selector, 
        all = rv_primers$evaluated_primers, 
        filtered = current.filtered.primers(), 
        optimized = optimal.primers())
    validate(need(data, "No data available for plotting."))
    validate(need(current.settings(), "Cannot plot coverage constraints since no settings were loaded."))
    # select only available constraints
    active.constraints <- input$selected_cvg_constraints
    #if (length(active.constraints) == 0) {  
        ## plot all constraints if none are selected to show something directly
        #active.constraints <- names(openPrimeR:::cvg_constraints(current.settings()))
    #}
    p <- openPrimeR:::plot_cvg_constraints(data, current.settings(), active.constraints)
    return(p)
})


output$constraint_plots_no_mismatches <- renderPlot({
    # plot number of mismatch binding events
    data <- switch(input$set_meta_selector, all = rv_primers$evaluated_primers, 
        filtered = current.filtered.primers(), optimized = optimal.primers())
    p <- openPrimeR:::plot_constraint.histogram.nbr.mismatches(data, allowed_nbr_of_mismatches())
    return(p)
})

eval.plot.width <- reactive({
    # width for evaluation plot
    primers <- switch(input$set_meta_selector, 
                      all = rv_primers$evaluated_primers, 
                      filtered = current.filtered.primers(), 
                      optimized = optimal.primers())
    if (length(primers) == 0) {
        return(1200)
    }
    idx <- grep("^EVAL_", colnames(primers))
    nbr.con <- length(idx)
    width <- openPrimeR:::get.plot.height(nbr.con, 40, 300)
    return(width)
})

eval.plot.height <- reactive({
    # height for evaluation plot
    primers <- switch(input$set_meta_selector, 
                      all = rv_primers$evaluated_primers, 
                      filtered = current.filtered.primers(), 
                      optimized = optimal.primers())
    if (length(primers) == 0) {
        return(900)
    }
    nbr.primers <- nrow(primers)
    height <- openPrimeR:::get.plot.height(nbr.primers, 25, 300)
    return(height)
})
output$constraint_fulfillment_plot <- renderPlot({
    primers <- switch(input$set_meta_selector, 
                      all = rv_primers$evaluated_primers, 
                      filtered = current.filtered.primers(), 
                      optimized = optimal.primers())
    validate(need(primers, "Please evaluate the primers first."))
    validate(need(nrow(primers) != 0, "No primers available."))
    input.primers <- list(primers)
    names(input.primers) <- primers$Run[1]
    validate(need(current.seqs(), "Please upload a set of templates first."))
    p <- openPrimeR:::plot_constraint_fulfillment(input.primers, 
            current.settings(),
            plot.p.vals = TRUE)
    return(p)
})
output$constraint_deviations <- renderPlot({
    # plot of deviations from target constraints
    constraint.df <- switch(input$set_meta_selector, 
                      all = rv_primers$evaluated_primers, 
                      filtered = current.filtered.primers(), 
                      optimized = optimal.primers())
    validate(need(constraint.df, "Please evaluate the primers first."))
    openPrimeR:::plot_constraint_deviation(constraint.df, current.settings())
})

output$constraint_stats <- renderPlot({
    # plot overview of fulfilled/failed constraints
    constraint.df <- switch(input$set_meta_selector, 
                      all = rv_primers$evaluated_primers, 
                      filtered = current.filtered.primers(), 
                      optimized = optimal.primers())
    validate(need(constraint.df, "Please evaluate the primers first."))
    constraint.settings <- NULL
    openPrimeR:::plot_constraint_fulfillment(constraint.df, current.settings())
})#, width = eval.plot.width, height = eval.plot.height)


########
# Shiny server functionalities for the comparison of primer sets
#############

# rv_comparison.data: template_path: paths to template comparison sets
    # primer_path: paths to primer comparison sets primers and seqs: lists with each
    # template.df for each primer (if only one template.df -> try to use it for all primers),
    # constraints: constraint settings that are valid for each or a single primer set
    # necessary) comp_primers and comp_templates: currently active sets used to
    # display comparison plots
rv_comparison.data <- reactiveValues(template_path = NULL, 
    primer_path = NULL, 
    primers = NULL,  # all loaded primers 
    seqs = NULL,  # all loaded templates
    primers_filtered = NULL, # selected filtered primers
    seqs_filtered = NULL, # selected filtered templates 
    constraints = NULL, primer_fnames = NULL, 
    comp_primers = NULL, comp_templates = NULL) # selected all primers, selected all templates

eval_comparison_con_width <- reactive({
    # width for region comparison plot
    # nbr of facets times facet width
    openPrimeR:::get.plot.height(2, 400)
})
eval_comparison_con_height <- reactive({
    # height for region comparison plot

    openPrimeR:::get.plot.height(ceiling(length(input$selected_other_plot) / 2), 400)
})
eval_comparison_cvg_height <- reactive({
    # height for region comparison plot
    nbr <- ifelse(length(input$selected_cvg_comp_constraints) == 0, 2, length(input$selected_cvg_comp_constraints))

    openPrimeR:::get.plot.height(ceiling(nbr / 2), 400)
})

eval_comparison_width <- reactive({
    # width for region comparison plot
    # nbr of facets times facet width
    openPrimeR:::get.plot.height(3, 600)
})
eval_comparison_height <- reactive({
    # height for region comparison plot
    openPrimeR:::get.plot.height(ceiling(length(plot.comp.primers())) / 3, 250)
})
region_comparison_width <- reactive({
    # width for region comparison plot
    # nbr of facets times facet width
    openPrimeR:::get.plot.height(3, 200)
})
region_comparison_height <- reactive({
    # height for region comparison plot
    openPrimeR:::get.plot.height(length(plot.comp.primers()) / 3, 300)
})
output$cvg_vs_size_plot <- renderPlot({
    
    validate(need(plot.comp.primers(), "Please upload primer sets for comparison and click on the Compare button first."))
    validate(need(plot.comp.templates(), "No template sequences corresponding to the primers were uploaded."))
    validate(need(length(plot.comp.templates()) == length(plot.comp.primers()),
        "Number of uploaded template sets did not agree with the number of uploaded primer sets."))
    p <- openPrimeR::plot_cvg_vs_set_size(plot.comp.primers(), plot.comp.templates(), show.labels = TRUE)
    return(p) 
    # plotly isn't used anymore (not rendered in Docker container somehow ..)
    #plotly::ggplotly(p) %>% plotly::layout(dragmode = "select")
})
#output$hover <- renderPrint({
    #d <- plotly::event_data("plotly_hover")
    #if (is.null(d)) "Hover events appear here (unhover to clear)" else d
  #})
#
#output$click <- renderPrint({
    #d <- plotly::event_data("plotly_click")
    #if (is.null(d)) "Click events appear here (double-click to clear)" else d
#})
#
#output$brush <- renderPrint({
    ##d <- plotly::event_data("plotly_selected")
    #if (is.null(d)) "Click and drag events (i.e., select/lasso) appear here (double-click to clear)" else d
#})
#
#output$zoom <- renderPrint({
    #d <- plotly::event_data("plotly_relayout")
    #if (is.null(d)) "Relayout (i.e., zoom) events appear here" else d
#})
output$comparison_plot_cvg <- renderPlot({
    # comparison of primer coverage for different primer sets
    validate(need(plot.comp.primers(), "Please upload primer sets for comparison and click on the Compare button first."))
    validate(need(plot.comp.templates(), "No template sequences corresponding to the primers were uploaded."))
    validate(need(length(plot.comp.templates()) == length(plot.comp.primers()), 
        "Number of uploaded template sets did not agree with the number of uploaded primer sets."))
    openPrimeR:::plot_template_cvg(plot.comp.primers(), plot.comp.templates())
})
output$comparison_plot_regions <- renderPlot({
    # binding regions comparison plot
    validate(need(plot.comp.primers(), "Please upload primer sets for comparison and click on the Compare button first."))
    validate(need(plot.comp.templates(), "No template sequences corresponding to the primers were uploaded."))
    validate(need(length(plot.comp.templates()) == length(plot.comp.primers()), 
        "Number of uploaded template sets did not agree with the number of uploaded primer sets."))
    group <- "all"
    direction <- "both"
    relation <- input$primer_comparison_relation
    openPrimeR:::plot_primer_binding_regions(plot.comp.primers(), plot.comp.templates(),
        direction, group, relation)
}, width = region_comparison_width, height = region_comparison_height, units = "px")

comparison_primer_choices <- reactive({
    # The primer sets that are available in the tool, depending on the chosen locus
    # of the user (input$template_comparison_locus)
    if (length(input$template_comparison_locus) == 0 || input$template_comparison_locus == "") {
        return(NULL)
    }
    options <- comparison.primer.choices(input$template_comparison_locus)
    return(options)
})
output$comp_cvg_constraints <- renderPlot({
    # plot cvg constraints
    validate(need(plot.comp.primers(), "Please upload primer sets for comparison and click on the Compare button first."))
    # select only the required cvg constraints:
    active.constraints <- input$selected_cvg_comp_constraints
    #if (length(active.constraints) == 0) {
        #active.constraints <- names(openPrimeR:::cvg_constraints(current.settings()))
    #}
    p <- openPrimeR:::plot_cvg_constraints(plot.comp.primers(), current.settings(), active.constraints)
    return(p)
}, width = eval_comparison_con_width, height = eval_comparison_cvg_height, units = "px")

availablePrimerComparisonUpdater <- observeEvent(input$template_comparison_locus, 
    {
    # update the primer set selector for comparison when selected template locus
    # (input$template_comparison_locus) changes
    # reset reactive values:
    reset.reactive.values(rv_comparison.data)
    choices <- comparison_primer_choices()
    updateSelectInput(session, "selected_comparison_primers", choices = choices)
    # update analysis identifier
    updateTextInput(session, "sample_name", 
        value = input$template_comparison_locus) 
})

comparison.table <- reactive({
    # overview of all loaded templates and primer sets for comparison
    primers <- rv_comparison.data$primers  # display all loaded primers
    templates <- rv_comparison.data$seqs
    validate(need(length(primers) != 0 || length(templates) != 0, "No data was uploaded for comparison."))
    out <- openPrimeR:::build.comparison.table(primers, templates)
    return(out)
})

output$uploaded_comp_data <- DT::renderDataTable({
    # shows the table with all data (primers/templates) available for comparison
    DT::datatable(comparison.table(), caption = "Overview of loaded primer-template pair sets. To compare only a subset of loaded primer sets, just select the corresponding rows in the table and click on \"Analyze->Compare\".", 
        options = list(searching = TRUE, processing = FALSE))
})

ComparisonPrimerInputObserver <- observeEvent(input$comparison_file, {
    # load primer comparison sets uploaded by the user as csv
    rv_values$comparison_primer_path <- input$comparison_file
})

ComparisonPrimerSuppliedObserver <- observeEvent(input$selected_comparison_primers, 
    {
        # update the path to the comparison primers when user selects a supplied primer
        # set
        if (length(input$selected_comparison_primers) == 0 || input$selected_comparison_primers == 
            "") {
            return(NULL)
        }
        choice.table <- data.frame(datapath = input$selected_comparison_primers, 
            name = NA, stringsAsFactors = FALSE)
        rv_values$comparison_primer_path <- choice.table
    })

read.compare.primers <- reactive({
    # loads the currently selected comparison primer sets
    if (length(rv_values$comparison_primer_path) != 0) {
        data <- withWarnings(openPrimeR:::read_primers(rv_values$comparison_primer_path$datapath))
        for (i in seq_along(data$warnings)) {
            warning <- data$warnings[[i]]
            message(warning)
        }
        for (i in seq_along(data$errors)) {
            error <- data$errors[[i]]
            print(error)
            if (inherits(error, "TemplateFormatIncorrect")) {
                toggleModal(session, "TemplateFormatIncorrect")
            }
        }
        if (length(rv_values$comparison_primer_path$datapath) == 1) {
            # ensure that loaded template set is always a list
            nam <- ifelse(nrow(data$value) == 0, rv_values$comparison_primer_path$name, data$value$Run[1])
            data <- list(data$value)
            names(data) <- nam
        } else {
            data <- data$value
        }
        updateTabsetPanel(session, "main", selected = "Comparison")
    }
    return(data)
})

comparisonTemplateOtherObserver <- observeEvent(input$comparison_templates, {
    # updates the path to the comparison templates upon user upload of template csv
    rv_comparison.data$comparison_template_path <- input$comparison_templates
})
comparisonTemplateIMGTObserver <- observeEvent(input$template_comparison_locus, {
    # updates the path to the comparison templates when user selects a supplied
    # template set
    rv_comparison.data$comparison_template_path <- get.supplied.comparison.template.path(input$template_comparison_locus)
})

read.comparison.templates <- reactive({
    # load templates for comparison
    if (length(rv_comparison.data$comparison_template_path) == 0) {
        return(NULL)
    }
    data <- withWarnings(openPrimeR:::read_templates(rv_comparison.data$comparison_template_path$datapath))
    for (i in seq_along(data$warnings)) {
        warning <- data$warnings[[i]]
        message(warning)
    }
    for (i in seq_along(data$errors)) {
        error <- data$errors[[i]]
        print(error)
        if (inherits(error, "TemplateFormatIncorrect")) {
            toggleModal(session, "TemplateFormatIncorrect")
        }
    }
    if (length(rv_comparison.data$comparison_template_path$datapath) == 1) {
        # ensure that loaded template set is always a list
        nam <- ifelse(nrow(data$value) == 0, rv_comparison.data$comparison_template_path$name, data$value$Run[1])
        data <- list(data$value)
        names(data) <- nam
    } else {
        data <- data$value
    }
    updateTabsetPanel(session, "main", selected = "Comparison")
    updateTextInput(session, "sample_name",  # update analysis identifier
                    value = rv_comparison.data$comparison_template_path$name)
    #runs.available <- sapply(data, function(x) "Run" %in% colnames(x))
    #validate(need(all(runs.available), "Input csv data is not supported. Please use the raw downloaded csv data."))
    return(data)
})
comparisonFileObserver <- observeEvent(rv_values$comparison_primer_path, {
    # input of primer sets for comparison stores all primers in
    # rv_comparison.data$primer_fnames for loading upon clicking the comparison
    # button
    input.primers <- rv_values$comparison_primer_path
    if (length(input.primers) != 0) {
        data <- read.compare.primers()
        l <- length(rv_comparison.data$primers)
        ## store all of the uploaded files in one variable and only remove if it is resetted
        if (length(rv_comparison.data$primers) != 0) {
            rv_comparison.data$primers <- c(rv_comparison.data$primers, data)
        } else {
            rv_comparison.data$primers <- data
        }
        names(rv_comparison.data$primers[(l + 1):length(rv_comparison.data$primers)]) <- names(data)
        if (length(rv_comparison.data$primer_fnames) != 0) {
            rv_comparison.data$primer_fnames <- c(rv_comparison.data$primer_fnames, 
                input.primers$name)
        } else {
            rv_comparison.data$primer_fnames <- input.primers$name
        }
        # switch to data overview table
        updateTabsetPanel(session, "main", selected = "Comparison")
        updateTabsetPanel(session, "selected_comparison_plot", selected = "loaded_data")
    }
})
comparisonTemplateObserver <- observeEvent(rv_comparison.data$comparison_template_path, 
    {
        # input template sets for comparison stores all template sets in
        # rv_comparison.data$seqs
        input.templates <- rv_comparison.data$comparison_template_path
        if (length(input.templates) != 0) {
            data <- read.comparison.templates()
            ## store all of the uploaded files in one variable and only remove if it is
            ## resetted
            if (length(rv_comparison.data$seqs) != 0) {
                rv_comparison.data$seqs <- c(rv_comparison.data$seqs, data)  # append to list
            } else {
                rv_comparison.data$seqs <- data  # start new list
            }
        # switch to data overview table
        updateTabsetPanel(session, "main", selected = "Comparison")
        updateTabsetPanel(session, "selected_comparison_plot", selected = "loaded_data")

        }
    })
current.comp.primers <- reactive({
    primers <- rv_comparison.data$primers
    # adjust primers according to selected rows in table 
    # if no rows are selected -> use all
    comp.table <- comparison.table()
    idx <- NULL
    if (length(comp.table) != 0 && nrow(comp.table) != 0) {
        idx <- as.numeric(rownames(comp.table))[input$uploaded_comp_data_rows_selected]
    } 
    if (length(idx) == 0) {
        idx <- seq_along(primers)
    }
    if (length(primers) != 0) {
        primers <- primers[idx]
    }
    return(primers)
})
current.comp.seqs <- reactive({
    seqs <- rv_comparison.data$seqs
    primers <- rv_comparison.data$primers
    if (length(seqs) == 1) {
        # use one template.df for all primer.data
        seqs <- replicate(length(primers), seqs[[1]], simplify = FALSE)
    }
    # adjust  seqs according to selected rows in table 
    # if no rows are -> use all entries
    comp.table <- comparison.table()
    idx <- NULL
    if (length(comp.table) != 0 && nrow(comp.table) != 0) {
        idx <- as.numeric(rownames(comp.table))[input$uploaded_comp_data_rows_selected]
    } 
    if (length(idx) == 0) {
        idx <- seq_along(seqs)
    }
    if (length(primers) != 0) {
        seqs <- seqs[idx]
        primers <- primers[idx]
        # set template cvg according to the primer sets:
        seqs <- lapply(seq_along(seqs), function(x) {
                            openPrimeR::update_template_cvg(seqs[[x]], primers[[x]])
        })
    }
    return(seqs)
})

primerComparisonObserver <- observeEvent(c(input$compare_primers), {
    # upon selecting input$compare_primers, retrieve the selected primers and set
    # them in rv_comparison.data$comp_primers and rv_comparison.data$comp_templates
    #
    # (TODO: change this -> set_meta_selector should be used to select the active set, not here!!!) ?
    #
    # Ensure that the run identifiers are unique before starting the analyses: 'set.run.names' does this
    primers <- openPrimeR:::set.run.names(current.comp.primers())
    seqs <- openPrimeR:::set.run.names(current.comp.seqs())
    # update comparison selection
    isolate({switch.view.selection("all", input$main, session)})
    updateTabsetPanel(session, "main", selected = "Comparison")
    updateTabsetPanel(session, "selected_comparison_plot", selected = "coverage_overview")
    # set in reactiveValues
    rv_comparison.data$comp_primers <- primers
    rv_comparison.data$comp_templates <- seqs
})
output$comparison_stats <- DT::renderDataTable({
    # data table giving an overview of the constraint stats of the compared primer
    # sets
    validate(need(plot.comp.primers(), "Please upload primer sets for comparison and click on the Compare button first."))
    validate(need(plot.comp.templates(), "No template sequences corresponding to the primers were uploaded."))
    validate(need(length(plot.comp.templates()) == length(plot.comp.primers()), 
        "Number of uploaded template sets did not agree with the number of uploaded primer sets."))
    data <- openPrimeR:::get_cvg_stats(plot.comp.primers(), plot.comp.templates(), for.viewing = TRUE)
    DT::datatable(data, options = list(searching = FALSE, processing = FALSE), caption = "Coverage statistics of the loaded primer sets.")
})

output$comparison_plot_deviation <- renderPlot({
    validate(need(plot.comp.primers(), "Please upload primer sets for comparison and click on the Compare button first."))
    openPrimeR:::plot_constraint_deviation(plot.comp.primers(), current.settings())
    # TODO: width and height -> width depends on nbr of sets?
})
output$comparison_plot_box <- renderPlot({
    # comparison plot for coverage
    validate(need(plot.comp.primers(), "Please upload primer sets for comparison and click on the Compare button first."))
    validate(need(plot.comp.templates(), "No template sequences corresponding to the primers were uploaded."))
    validate(need(length(plot.comp.templates()) == length(plot.comp.primers()), 
        "Number of uploaded template sets did not agree with the number of uploaded primer sets."))
    openPrimeR:::plot.comparison.box(plot.comp.primers(), plot.comp.templates())
})
output$comparison_plot_constraint <- renderPlot({
    validate(need(plot.comp.primers(), "Please upload primer sets for comparison and click on the Compare button first."))
    p <- openPrimeR:::plot_constraint(plot.comp.primers(),
                        current.settings(), input$selected_other_plot)
    return(p)
}, width = eval_comparison_con_width, height = eval_comparison_con_height, units = "px")

output$comparison_plot_mismatches <- renderPlot({
    # comparison plot for the number of mismatches in primer sets
    validate(need(plot.comp.primers(), "Please upload primer sets for comparison and click on the Compare button first."))
    validate(need(plot.comp.templates(), "No template sequences corresponding to the primers were uploaded."))
    validate(need(length(plot.comp.templates()) == length(plot.comp.primers()), 
        "Number of uploaded template sets did not agree with the number of uploaded primer sets."))
    p <- openPrimeR:::plot_primer.comparison.mismatches(plot.comp.primers(), plot.comp.templates(),
        allowed_nbr_of_mismatches())
    return(p)
})
output$comparison_plot_evaluation_ui <- renderUI({
    # ui output of evaluation to prevent overplotting
    plotOutput("comparison_plot_evaluation",
        width = paste0(region_comparison_width(), "px"), 
        height = paste0(region_comparison_height(), "px"))
})

output$comparison_plot_evaluation <- renderPlot({
    # comparison plot regarding constraint evaluation evaluate primers according to
    # active constraints
    validate(need(plot.comp.primers(), "Please upload primer sets for comparison and click on the Compare button first."))
    # add melting_temp_diff manually here
    # TODO: signature plot disabled until ggiraph bug fixed in the pkg
    p <- openPrimeR:::plot_constraint_fulfillment(plot.comp.primers(),
            current.settings(), plot.p.vals = FALSE)
    #p <- openPrimeR:::plot_constraint_signature(plot.comp.primers(),
            #plot.comp.templates(),
            #constraint.settings, plot.p.vals = FALSE, ncol = NULL)
    # note: use active settings here to make the plot dynamic
    return(p)
}, width = region_comparison_width, height = region_comparison_height, units = "px")
#  fix width for radar plot: , width = 1200, height = 1200, units = "px")

output$comparison_overview_table <- DT::renderDataTable({
    # overview of constraints in comparison primer sets
    validate(need(plot.comp.primers(), "Please upload primer sets for comparison and click on the Compare button first."))
    validate(need(plot.comp.templates(), "No template sequences corresponding to the primers were uploaded."))
    validate(need(length(plot.comp.templates()) == length(plot.comp.primers()), 
        "Number of uploaded template sets did not agree with the number of uploaded primer sets."))
    table <- openPrimeR:::get_comparison_table(plot.comp.templates(), plot.comp.primers())
    return(DT::datatable(table, 
        caption = paste("Overview of compared primer sets.",
                  "Values in angular brackets indicate inter-quartile ranges."),
        rownames = FALSE,  extensions="Responsive"))
})

observeEvent(input$load_all_comparison_sets, {
    # load all available comparison primer sets reset current table and then load all
    rv_comparison.data$primers <- NULL
    rv_comparison.data$primers_filtered <- NULL
    rv_comparison.data$constraints <- NULL
    # load all primer set options
    choices <- comparison_primer_choices()
    choice.table <- data.frame(datapath = choices, name = names(choices), stringsAsFactors = FALSE)
    rv_values$comparison_primer_path <- choice.table
})

comparisonFilterObserver <- observeEvent(input$quickFilterButton_compare, {
    # filter comparison primer sets using the selected constraints Create a callback
    # try to set filtered so comparison compares the correct primers
    updateSelectInput(session, "set_meta_selector", selected = "filtered")
    # function to update progress.
    progress <- shiny::Progress$new()
    progress$set(message = "Evaluating", value = 0)
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
    ### 
    primers <- current.comp.primers()
    templates <- current.comp.seqs()
    if (length(templates) == 0 || length(primers) == 0) {
        return
    } 
    con <- constraints()
    active.constraints <- names(openPrimeR:::constraints(current.settings()))
    filter.data <- openPrimeR:::filter.comparison.primers(primers, templates, 
        active.constraints, current.settings(),
        updateProgress)
    # save data in reactive rv_values lists
    rv_comparison.data$primers_filtered <- filter.data$primers
    rv_comparison.data$seqs_filtered <- filter.data$templates
    
})
plot.comp.primers <- reactive({
    primers <- switch(input$set_meta_selector,
            "all" = rv_comparison.data$comp_primers,
            "filtered" = rv_comparison.data$primers_filtered,
            "optimized" = NULL) 
    #print("No primer sets:")
    #print(length(primers))
    return(primers)
})
plot.comp.templates <- reactive({
    templates <- switch(input$set_meta_selector,
            "all" =  rv_comparison.data$comp_templates,
            "filtered" =  rv_comparison.data$seqs_filtered,
            "optimized" = NULL)
    #print("no template sets:")
    #print(length(templates))
    return(templates)
})
comparisonResetObserver <- observeEvent(input$reset_rv_comparison.data, {
    # reset the rv_comparison.data (loaded primers, templates, constraints)
    reset.reactive.values(rv_comparison.data)
    session$sendCustomMessage(type = "resetFileInputHandler", "comparison_file")
    session$sendCustomMessage(type = "resetFileInputHandler", "comparison_templates")
    session$sendCustomMessage(type = "resetFileInputHandler", "comparison_constraint_files")
})

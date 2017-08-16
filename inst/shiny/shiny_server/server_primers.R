############
# Shiny server functionalities relating to the input/evaluation of primers
############
# REACTIVE VALUES
# rv_primers: primer data frames at different analysis stages

rv_primers <- reactiveValues("PrimerTab" = NULL, # loaded primers
                             "PrimerTabFiltered" = NULL, # displayed filtered primers
                             "PrimerTabOptimized" = NULL, # displayed opti primers
                             "last_filter_type" = NULL, # last filtering type performed: either filtered or opti
                             "filtered" = NULL, # filtered primer data (list with stats etc)
                             "filtered_opti" = NULL, # filtered opti primer data
                             "evaluated_primers" = NULL, # evaluated primer set
                             "optimal_data" = NULL, # all data from optimization

                             "optimized" = NULL, # optimized primer set
                             "selected_idx" = NULL, # selected row in primer table. not used at the moment
                             "available_constraints" = NULL, # names of computed constraints by evaluation
                             "all" = NULL) # input data

input.primers <- reactive({
    # Loads the input primer data
    primerFile <- rv_cur.input.data$primers
    #print(paste0("input.primers(): primerFile is: ", primerFile))
    if (length(primerFile) == 0) {
        return(NULL)
    }
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Reading primers", value = 0)
    on.exit(progress$close())
    updateProgress <- function(value = NULL, detail = NULL, option = NULL) {
        if (is.null(value)) {
            value <- progress$getValue()
            value <- value + (progress$getMax() - value) / 5
        }
        if (option == "inc") {
            progress$inc(amount = value, detail = detail) 
        } else {
            progress$set(value = value, detail = detail)
        }
    }
    primers <- isolate({withWarnings(openPrimeR:::read_primers(
                primerFile$datapath, input$fw_primer_id, 
                input$rev_primer_id, input$use_ambig, 
                input$max_degeneracy, 
                updateProgress = updateProgress
                ))})
    for (i in seq_along(primers$errors)) {
        error <- primers$errors[[i]]
        print(error)
        if (inherits(error, "NotifyPrimersNoDirection")) {
            toggleModal(session, "NotifyPrimersNoDirection")
        } else if (inherits(error, "FastaAlphabetError")) {
            toggleModal(session, "FastaAlphabetError")   
        } else if (inherits(error, "TemplateFormatIncorrect")) {
            toggleModal(session, "TemplateFormatIncorrect")
        } else {
            toggleModal(session, "TemplateFormatIncorrect")
        }
    }
    for (i in seq_along(primers$warnings)) {
        warning <- primers$warnings[[i]]
        if (inherits(warning, "NotifyPrimersMissingKeyword")) {
            toggleModal(session, "NotifyPrimersMissingKeyword")
        } else if (inherits(warning, "NotifyPrimersDuplicateDirections")) {
            toggleModal(session, "NotifyPRimersDuplicateDirections")
        }
    }
    if (length(primers$errors) != 0) {
        session$sendCustomMessage(type = "resetFileInputHandler", 'primers_file')
        primers <- NULL
    } else {
        isolate({updateTextInput(session, "sample_name", value = paste0(primers$value$Run[1], "|", input$sample_name))}) # update analysis identifier w/ primer ID
        primers <- primers$value
    }
    if (length(primers) == 0) {
        rv_primers$PrimerTab <- NULL
        validate(need(primers, "No primers are available."), errorClass = "fatal")
    } 
    # select primer view in UI:
    if (isolate(!input$load_eval_primers || length(rv_primers$evaluated_primers) == 0)) {
        # only switch to the primers tab for raw primers or
        # when no primers haven't been loaded yet, otherwise stay
        # on current tab to quickly gauge different sets
        updateTabsetPanel(session, "main", selected = "Primers")
    }
    isolate({switch.view.selection("all", input$main, session)})
    return(primers)
})

output$primer_restriction_sites <- DT::renderDataTable({ 
    DT::datatable(primer.restriction.sites())
})
AdapterModalObserver <- observeEvent(input$check_adapters, {
    sites <- primer.restriction.sites() # compute primer restriction sites
    if (length(sites) != 0 && nrow(sites) != 0) {
        ## trigger modal for adapter display
            # modal for adapter check:
        toggleModal(session, "AdapterModal")  # show modal when selection changes
    } else {
        toggleModal(session, "NoAdapterModal")
    }
})
primer.restriction.sites <- reactive({
    template.df <- current.seqs()
    site.df <- NULL
    if (length(template.df) != 0) {
        # check for restriction sites
        # Create a Progress object
        progress <- shiny::Progress$new()
        progress$set(message = "Checking adapters", value = 0)
        on.exit(progress$close())
        updateProgress <- function(value = NULL, detail = NULL, option = NULL) {
            if (is.null(value)) {
                # if value is NULL, increase bar by 1/5th of the remaining distance to cover
                value <- progress$getValue()
                value <- value + (progress$getMax() - value) / 5
            }
            if (option == "inc") {
                progress$inc(amount = value, detail = detail) 
            } else {
                progress$set(value = value, detail = detail)
            }
        }
        site.df <- openPrimeR:::check_restriction_sites(input.primers(), template.df, 
                                adapter.action = "warn",
                                updateProgress = updateProgress)
    } 
    return(site.df)
})
primers.IMGT.view <- reactive({ 
    # Returns the names/paths of primer sets that are available when IMGT templates have been selected 
    if (length(input$IMGT_DB_locus) == 0) {
        return(NULL)
    }
    if (input$load_eval_primers) {
        # load pre-evaluated comparison sets
        primer.folder <- system.file("extdata", "IMGT_data", "comparison", "primer_sets", package=  "openPrimeR")
    } else {
        # load raw fasta seqs
        primer.folder <- system.file("extdata", "IMGT_data", "primers", package=  "openPrimeR")
    }
    available.primer.loci <- list.dirs(primer.folder, recursive = FALSE)
    idx <- grep(input$IMGT_DB_locus, available.primer.loci)
    if (length(idx) == 0) {
        return(NULL) # no primers available for current locus selection
    }
    path <- available.primer.loci[idx]
    files <- list.files(path, full.names = TRUE)
    if (length(files) == 0) {
        return(NULL)
    }
    primer.paths <- primer.set.choices(files)
    return(primer.paths)
})
selected.IMGT.primers <- reactive({
    # Returns the paths to the currently selcted IMGT primers
    sel.primers <- input$IMGT_primers
    if (sel.primers == "") { # nothing selected
        return(NULL)
    }
    out <- list("datapath" = sel.primers, "name" = "IMGT_primers")
    return(out)
})
IMGT_EvaluatedObserver <- observeEvent(c(input$load_eval_primers, input$primer_upload_choice), {
    # disable the iupac ambiguity action selector if we load evaluated primer sets
    if (input$load_eval_primers) {
        # loading processed csv -> can't disambiguate primers
        shinyjs::disable("use_ambig")
    } else {
        # loading raw fasta -> can disambiguate primers
        shinyjs::enable("use_ambig")
    }
})

availablePrimerUpdater <- observeEvent(c(input$IMGT_DB_locus, input$load_eval_primers), {
    # Updates the input$IMGT_primers field when the input locus changes based on the available primers determined by primers.IMGT.view()
    primer.choices <- primers.IMGT.view()
    if (is.null(primer.choices)) {
        primer.choices <- character(0) # remove all choices
    }
    updateSelectInput(session, "IMGT_primers", choices = primer.choices)
})

primer.data <- reactive({
    # Loads the input primers from an input fasta file.
    primers <- input.primers()
    if (length(primers) == 0) {
        return(NULL)
    }
    primer.options <- primers$ID
    primer.options <- primer.options[order(primer.options)]
    # update the list of available primers in the Coverage tab
    updateSelectInput(session, "selected_primer", choices = primer.options)
    return(primers)
})
primer_subset_out <- reactive({
    # only used when downloading a set
    # ensure that subset data frame has all properties annotated correctly.
    # (some properties are not refreshed in the subset function)
    if (input$selected_subset_size == "") { # no subset selected
        return(NULL)
    }
    subset.df <- primer_subset()
    con <- constraints()
    active.constraints <- con[["active"]]
    seqs <- current.seqs()
    settings <- current.settings() 
    # compute all properties again (annotation of binding positions etc.)
    subset.df <- openPrimeR:::check_constraints(subset.df, seqs,
                                settings, active.constraints, 
                                for.shiny = TRUE)
    return(subset.df)
})
primer_subset <- eventReactive(c(input$selected_subset_size, primer_subsets()), {
    # Internal function for primer subsets
    if (input$selected_subset_size == "" || length(primer_subsets()) == 0) {
        return(NULL)
    }
    subset <- primer_subsets()[[as.numeric(input$selected_subset_size)]]
    return(subset)
}, ignoreNULL = FALSE)

output$primer_subset_table <- DT::renderDataTable({ 
    # Returns a table for the primer subset of the selected size 
    validate(need(primer_subset(), "No subset computed yet."))
    df <- view.subset.primers(primer_subset(), current.seqs(), run.mode(), input$view_cvg_individual)
    DT::datatable(df, caption = "Primers of the optimal subset.", escape=FALSE, options = list(processing = FALSE))
})

primer_subsets <- reactive({
    # Create primer subsets for the selected primers in the coverage tab.
    k <- 1
    data <- switch(input$set_meta_selector,
            "all" = rv_primers$evaluated_primers,
            "filtered" = current.filtered.primers(),
            "optimized" = optimal.primers())
    validate(need(data, "Please compute the primer coverage first."))
    template.data <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    validate(need(template.data, "No templates available yet."))
    if (length(input$selected_group_coverage) == 0 || "all" %in% input$selected_group_coverage) {
        groups <- NULL
    } else {
        groups <- input$selected_group_coverage
    }
    subsets <- openPrimeR:::subset_primer_set(data, template.data, k, groups, NULL, NULL)
    # update subset slider:
    if (length(subsets) != 0) {
        # update subset selector: show cvg of each subset size
        cvg.string <- sapply(subsets, function(x) paste("Coverage ", round(openPrimeR::get_cvg_ratio(x, template.data), 2) * 100, "%", sep = ""))
        set_size <- nrow(subsets[[length(subsets)]]) # paired set size
        labels <- paste("Size ", seq_len(set_size), " (", cvg.string, ")", sep = "")
        opts <- seq_len(set_size)
        names(opts) <- labels
        updateSelectInput(session, "selected_subset_size", choices = opts)
    }
    return(subsets)
})

output$primer_subset_coverage <- renderPlot({
# Plots the coverage achieved by each primer subset
    subsets <- primer_subsets()
    template.data <- switch(input$set_meta_selector,
            "all" = rv_templates$cvg_all,
            "filtered" = rv_templates$cvg_filtered,
            "optimized" = rv_templates$cvg_optimized)
    openPrimeR:::plot_primer_subsets(subsets, template.data, 
                                    required.cvg = input$required_opti_cvg)
})

PrimerTabObserver <- observe({
    # sets the current PrimerTab and rv_primers rv_values. done in an observer such that view options are respected.
    ####
    # ALL DATA TAB:
    ####
    # n.b.: observer can't output validate result ..
    if (length(primer.data()) != 0) { # input primers
        rv_primers$PrimerTab <- openPrimeR:::view.input.primers(primer.data(), run.mode())
        #print("setting new primer tab according to primer.data()")
        rv_primers$all <- primer.data()
    } else { # 
        rv_primers$PrimerTab <- NULL
        rv_primers$all <- NULL
    }
    if (length(rv_primers$evaluated_primers) != 0 && length(current.seqs()) != 0 && length(run.mode()) != 0) { # evaluated primers
        rv_primers$PrimerTab <- view.filtered.primers.all(rv_primers$evaluated_primers, current.seqs(), run.mode(), input$view_cvg_individual)
    } 
    #####
    # FILTERED DATA TAB
    ######
    if (length(current.filtered.primers()) != 0 && length(current.seqs()) != 0 && length(run.mode()) != 0) {
        rv_primers$PrimerTabFiltered <- view.filtered.primers(current.filtered.primers(), current.seqs(), run.mode(), input$view_cvg_individual)
    }
    #####
    # OPTIMIZED DATA TAB
    #####
    if (length(optimal.primers()) != 0 && length(current.seqs()) != 0 && length(run.mode()) != 0) {
        opti <- optimal.primers()
        rv_primers$PrimerTabOptimized <- view.optimized.primers(opti, current.seqs(), run.mode(), input$view_cvg_individual)
    }
}, priority = 5) # set high priority for updates ..



optimal.primers <- reactive({
    # function for optimal primer data frame
    primer.data <- rv_primers$optimal_data
    if (length(primer.data) == 0) {
        return(NULL)
    }
    return(primer.data$opti)
})
cur_primer_detail <- reactive({ # think of a table layout? TODO or do something else?
    primers <- switch(input$set_meta_selector,
                "all" = rv_primers$PrimerTab, 
                "filtered" = rv_primers$PrimerTabFiltered,
                "optimized" = rv_primers$PrimerTabOptimized)
    if (length(rv_primers$selected_idx) == 0 || length(primers) == 0) {
        data <- NULL
    } else {
       data <- primers[rv_primers$selected_idx,]  
       # TODO: use a function to convert single primer values to a nice table
    }
    return(data)
})
##############
###############
#output$primer_detail_table <- DT::renderDataTable({
    ## select active primer table:
    #primers <- switch(input$set_meta_selector,
                #"all" = rv_primers$PrimerTab, 
                #"filtered" = rv_primers$PrimerTabFiltered,
                #"optimized" = rv_primers$PrimerTabOptimized)
    ## show a table with the properties of the currently selected primer
    #validate(need(primers, "No primers available."))
    #validate(need(cur_primer_detail(), "No details to show at the moment."))
    #DT::datatable(cur_primer_detail())
    ##DT::datatable(cur_dimer_detail(), caption = "", escape = FALSE, rownames = FALSE) %>% 
        ##DT::formatStyle("DeltaG", backgroundColor = styleInterval(cutoff, c("#ff9999", 
            ##"#99d6ff")), )
#})

#primerDetailObserver <- observeEvent(input$PrimerTab_rows_selected, { 
    # if primer is selected, show properties of primer
    # TODO?
    #if (FALSE) { # TODO: think about this feature!
        #sel.ID <- input$PrimerTab_row_last_clicked  # last clicked row: since ID is the first column, we need to match to ID
        ## store in reactive rv_values to access by reactive function
        #rv_primers$selected_idx <- as.numeric(sel.ID)  # only works if rownames are reset to 1:N
        #toggleModal(session, "PrimerDetail")  # show modal when selection changes
    #}
#})
###################
output$runModeText <- renderUI({
    text <- paste0("Coverage mode: ", run.mode())
    # Change to blue, size is ok?
    return(HTML(text))
})

output$designText <- renderUI({
    settings <- current.settings()
    allowed.mismatches <- openPrimeR::conOptions(settings)$allowed_mismatches
    run.mode <- input$design_direction
    init.mode <- input$init_algo
    opti.algo <- input$optimization_algorithm
    template.df <- current.seqs()    
    required.cvg <- input$required_opti_cvg
    text <- paste0(create.design.string(allowed.mismatches, run.mode, init.mode, opti.algo, template.df, required.cvg))
    # add a warning about the runtime
    text <- paste0(text, paste0(" Dependent on your data set, the computations may take a considerable amount of time (e.g. multiple hours).",
                    " The computations can only be interrupted by forcefully stopping the tool. Before designing a primer set, you may want to estimate whether it is possible to find a reasonable set of primers for the provided templates by evaluating the problem's difficulty."))
    return(HTML(text))
})
output$designTextDiff <- renderUI({
    problem.text <- ""
    if (length(problem.difficulty()) != 0) {
        a <- switch(problem.difficulty()$Classification,
                    "very_easy" = "a very easy",
                    "easy" = "an easy",
                    "medium" = "a typical",
                    "hard" = "a hard",
                    "very_hard" = "a very hard")
        interpretation <- NA
        if (grepl("easy", problem.difficulty()$Classification)) {
            interpretation <- paste("It should be possible to design",
                                "a small set of primers covering the template sequences.")
        } else if (grepl("medium", problem.difficulty()$Classification)) {
            interpretation <- paste("It should be possible to design",
                                    "a primer set, but it may be hard",
                                    "to cover all templates or reach the",
                                    "target coverage.")
        } else if(grepl("hard", problem.difficulty()$Classification)) {
            interpretation <- paste("It may be difficult to design",
                                    "a small set of primers covering",
                                    "the template sequences.",
                                    "Is it possible that there is a more",
                                    "conserved binding region you could choose?")
        } 
        if (problem.difficulty()$Uncertain) {
            interpretation <- paste(interpretation, 
                                    "Note that the fit of the beta distribution was", 
                                   "not good enough to allow for a ",
                                   "confident classification of the problem's difficulty.")
        }
        conf <- paste0(round(problem.difficulty()$Confidence * 100, 2), "%")
        diff.text <- paste0("This seems to be ", a, 
                           " primer design problem (confidence ",
                           conf, ").")
        problem.text <- paste(diff.text, interpretation)
        b.fw <- problem.difficulty()$Nbr_primers_fw
        b.rev <- problem.difficulty()$Nbr_primers_rev
        if (!is.na(b.fw) || !is.na(b.rev)) {
            if (length(input$design_direction) == 0) {
                b <- "unknown" 
            } else if (input$design_direction == "both") {
                b <- paste0(b.fw + b.rev, " primers (", b.fw, " fw / ", 
                            b.rev, " rev)")
            } else if (input$design_direction == "fw") {
                b <- paste0(b.fw, " primers")
            } else {
                b <- paste0(b.rev, " primers")
            }
            primer.text <- paste0("The estimated required number of primers is ", b, ".")
            problem.text <- paste(problem.text, primer.text)
       }
    }
    return(HTML(problem.text))
})

problem.difficulty <- eventReactive(input$evaluate_difficulty, {
    if (input$evaluate_difficulty == 0) {
        return(NULL)
    }
    design.diff <- withWarnings(openPrimeR::classify_design_problem(current.seqs(), 
                        input$design_direction, 
                        min(input$allowed_primer_length),
                        input$evaluate_difficulty_primers,
                        input$required_opti_cvg 
                        ))
    # check for warnings:
    for (i in seq_along(design.diff$warnings)) {
        warning <- design.diff$warnings[[i]]
        if (inherits(warning, "ProblemEstimationProblem")) {
            toggleModal(session, "ProblemEstimationProblem")
        }
        warning(warning)
    }
    design.diff <- design.diff$value
    if (length(design.diff) == 0) {
        # no result since distribution couldn't be estimated.
        return(NULL)
    }
    # change traffic light
    active.class <- NA
    if (grepl("easy", design.diff$Classification)) {
        active.class <- "green"
    } else if (grepl("medium", design.diff$Classification)) {
        active.class <- "orange"
    } else if (grepl("hard", design.diff$Classification)) {
        active.class <- "red"
    } else {
        active.class <- NULL
    }
    if (design.diff$Uncertain) {
        # don't show a light when uncertain.
        active.class <- NULL
    }
    # 1. set/remove active classes
    classes <- c("green", "orange", "red")
    for (i in seq_along(classes)) {
        selector <- paste0("#light #", classes[i])
        if (classes[i] %in% active.class) {
            shinyjs::addCssClass(class = "active", selector = selector)
        } else {
            shinyjs::removeCssClass(class = "active", selector = selector)
            
        }
    }
    shinyjs::show(selector = "#light")
    return(design.diff)
}, ignoreNULL = FALSE)

optimal.primer.data <- observeEvent(input$optimizeButton, {
    if (input$optimizeButton == 0) {
        return(NULL)
    }
    # close design verification pop-up when "go" button is pressed
    toggleModal(session, "DesignVerification")
    # create a modal if no data is available
    if (length(current.seqs()) == 0) {
        session$sendCustomMessage(type='jsCode', list(value = "$('#NotifyNoDataAvailable').modal('show')"))
        return(NULL)
    }
    # Create a Progress object
    progress <- shiny::Progress$new()
    progress$set(message = "Optimizing", value = 0)
    on.exit(progress$close())
    # Create a callback function to update progress.
    updateProgress <- function(value = NULL, detail = NULL, option = NULL) {
        if (is.null(value)) {
            # if value is NULL, increase bar by 1/5th of the remaining distance to cover
            value <- progress$getValue()
            value <- value + (progress$getMax() - value) / 5
        }
        if (option == "inc") {
            progress$inc(amount = value, detail = detail) 
        } else {
            progress$set(value = value, detail = detail)
        }
    }
    constraint.settings <- constraints()$active_settings
    required.cvg <- input$required_opti_cvg
    if (!input$relaxation_active) {
        # deactivate the relaxation procecdure by setting target cvg to 0
        required.cvg <- 0
    }
    cur.results.loc <- NULL # change to directory for debugging of results
    settings <- current.settings()
    primer.data <- withWarnings(openPrimeR:::design_primers(current.seqs(), 
                    input$design_direction, settings, 
                    input$init_algo, input$optimization_algorithm, 
                    required.cvg = input$required_opti_cvg, 
                    timeout = Inf,
                    max.degen = input$max_degeneracy, 
                    conservation = input$required_conservation, 
                    sample.name = input$sample_name, 
                    cur.results.loc = cur.results.loc,
                    updateProgress = updateProgress))
    # check for warnings:
    for (i in seq_along(primer.data$warnings)) {
        warning <- primer.data$warnings[[i]]
        if (inherits(warning, "AllowedRegionTooShort")) {
            toggleModal(session, "AllowedRegionTooShort")
        }
        warning(warning)
    }
    # check for errors:
    for (i in seq_along(primer.data$errors)) {
        error <- primer.data$errors[[i]]
        if (inherits(error, "PrimersDuplicateDirections")) {
            toggleModal(session, "NotifyPrimersDuplicateDirections")
        }
        print(error)
    }
    primer.data <- primer.data$value
    rv_primers$optimized <- primer.data$opti
    rv_primers$optimal_data <- primer.data # the optimal primer data
    rv_primers$filtered_opti <- primer.data$filtered # filtering data (not only primers, but also stats)
    filtered.templates <- openPrimeR:::update_template_cvg(current.seqs(), primer.data$filtered$data, run.mode()) # update templates with cvg info for filtered primers
    rv_templates$cvg_filtered <- filtered.templates
    rv_primers$last_filter_type <- "opti"
    # determine whether constraints were relaxed
    any.relaxed <- rep(FALSE, 2)
    for (i in seq_along(primer.data$used_constraints)) {
        # for fw/rev constraints:
        used.constraints <- openPrimeR::constraints(primer.data$used_constraints[[i]])
        opti.relaxed <- openPrimeR:::were.constraints.relaxed(used.constraints, openPrimeR::constraints(settings))
        any.relaxed[i] <- opti.relaxed
    }
    if (any(any.relaxed)) {
        rv_values$relax_info <- "$('#RelaxInfoOpti').modal('show')"
    } else {
        rv_values$relax_info <- NULL
    }
    rv_templates$cvg_optimized <- openPrimeR:::update_template_cvg(current.seqs(), primer.data$opti, run.mode()) # update templates with cvg info
    switch.view.selection("optimized", input$main, session) # switch to optimized primer view
})

output$PrimerTab <- DT::renderDataTable({ 
# Output a table showing the currently selected primers according to input$set_meta_selector
    # render the currently selected primer table
    withProgress(message = 'Rendering primer table ...', value = 0, {

        data <- switch(input$set_meta_selector,
                "all" = rv_primers$PrimerTab, # current primer table, but no exclusion
                "filtered" = rv_primers$PrimerTabFiltered,
                "optimized" = rv_primers$PrimerTabOptimized
                )
        }
    )
    #print("Rendered primer tab is:" )
    #print(summary(data))
    validate(need(data, "There aren't any available primers. Please check your input files and settings.")) 
    DT::formatStyle(DT::datatable(asS3(data), caption = "Overview of the primers.", 
                        escape=FALSE, options = list(processing = FALSE), 
                        extensions = "Responsive"),
                    'Direction', backgroundColor = DT::styleEqual(c("fw", "rev", "both"), c('#f4f6f7', '#f7fffa', '#edfffe')), target = "row")
})

notifyRelaxation <- observe({
    # show bsmodal with info about relaxed constraints during the optimization to the user. Depends on rv_values$relax_info being set after the optimization.
    if (length(rv_values$relax_info) != 0) {
        session$sendCustomMessage(type='jsCode', list(value = rv_values$relax_info))
    }
})
notifyNotAllowedBinding <- observe({
    # notification when the number of primers binding in disallowed regions exceeds the allowed ratio after primer evaluation.
    primers <- rv_primers$evaluated_primers
    if (length(primers) != 0 && "primer_coverage" %in% colnames(primers)) {
        check.allowed.binding <- strsplit(primers$Binding_Region_Allowed, split = ",")
        disallowed.binding.primer.idx <- which(sapply(check.allowed.binding, function(x) "Disallowed" %in% x))
        nbr.found <- length(disallowed.binding.primer.idx)
        nbr.allowed <- input$allowed_other_binding_ratio * nrow(primers)
        if (nbr.found > nbr.allowed) {
            session$sendCustomMessage(type='jsCode', list(value = "$('#NotifyNotAllowedBinding').modal('show')"))
        }
    }
})
InputPrimerObserver <- observeEvent(input$primer_file, { 
    # update current input primer file on user upload of primers
    #print(paste("Uploaded primer file: ", innput$primer_file))
    rv_cur.input.data$primers <- input$primer_file
})
IMGT_PrimerObserver <- observeEvent(input$IMGT_primers, {
    # update current input primer file on user selection of provided IMGT primers
    if (input$IMGT_primers == "") {
        # Dont update here on empty selection 
        return()
    }  
    #print(paste("IMGT primer file: ", selected.IMGT.primers()))
    rv_cur.input.data$primers <- selected.IMGT.primers()
})
run.mode <- reactive({
    # get analysis mode for primers. either fw/rev/both, depending on the directionality of the primers.
    run.mode <- openPrimeR:::get.analysis.mode(primer.data())
    if (input$primer_analysis_type == "design" || is.null(run.mode)) {
        run.mode <- input$design_direction
    }
    return(run.mode)
})
observeEvent(input$template_scenario, {
    # load evaluated primers if available templates are selected
    # otherwise, load non-evaluated primer sets
    if (input$template_scenario == "supplied") {
        updateCheckboxInput(session, "load_eval_primers", value = TRUE)
    } else {
        updateCheckboxInput(session, "load_eval_primers", value = FALSE)
    }
})

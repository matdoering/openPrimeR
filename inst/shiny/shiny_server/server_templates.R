######################
# Shiny server functionalities for the input of templates and allowed regions
######################

# REACTIVE VALUES
rv_templates <- reactiveValues("SeqTab" = NULL, # displayed sequences
                              "SeqTabFiltered" = NULL, # updated cvg from filtering/optimization
                              "SeqTabOptimized" = NULL, # templates with cvg from optimization primer set
                              "cvg_all" = NULL, # annotated template cvg for all primers
                              "cvg_filtered" = NULL, # annotated template cvg for filtered primers
                              "cvg_optimized" = NULL, # annotated template cvg for optimized primers
                              "raw_seqs" = NULL, # all seqs without any changes (no view)
                              "load_IMGT_templates" = FALSE, # boolean to indicate whether IMGT templates are to be loaded.
                              "cur_cvg_regions" = NULL # the allowed binding regions for which cvg has been computed for 'all'primers
)

InputDataObserverTemplates <- observeEvent(input$sequence_file, {
    # Updates rv_cur.input.data reactive values when templates are uploaded manually.
    #
    #   input$sequence_file: the template fasta input
    if (length(rv_cur.input.data$templates_exon) != 0 && length(rv_cur.input.data$templates_leader) != 0) {
        # reset previous input data
        rv_cur.input.data$templates_leader <- NULL
    }
    rv_templates$load_IMGT_templates <- FALSE # no IMGT input
    rv_cur.input.data$templates_exon <- input$sequence_file
})
InputDataObserverLeaders_fw <- observeEvent(input$leader_file, { 
    # Updates rv_cur.input.data reactive values when fw individual binding regions are uploaded manually.
    #
    #   input$leader_file: the individual binding region file   
    if (length(rv_cur.input.data$templates_exon) != 0 && length(rv_cur.input.data$templates_leader) != 0) {
        # reset previous input data
        rv_cur.input.data$templates_exon <- NULL
    }
    rv_cur.input.data$templates_leader <- input$leader_file
})
InputDataObserverLeaders_rev <- observeEvent(input$leader_file_rev, { 
    # Updates rv_cur.input.data reactive values when rev individual binding regions are uploaded manually.
    #
    #   input$leader_file: the individual binding region file   
    if (length(rv_cur.input.data$templates_exon) != 0 && length(rv_cur.input.data$templates_leader_rev) != 0) {
        # reset previous input data
        rv_cur.input.data$templates_exon <- NULL
    }
    rv_cur.input.data$templates_leader_rev <- input$leader_file_rev
})
seq.data.input <- reactive({
    # Loads template data using "rv_cur.input.data" reactiveValues file.
    seqFile <- rv_cur.input.data$templates_exon
    if (is.null(seqFile)) {
      return(NULL)
    }
    is.IMGT.data <- isolate(rv_templates$load_IMGT_templates)
    rm.keywords <- NULL
    rm.duplicated <- FALSE
    if (is.IMGT.data) {
        # overwrite header structure and id column for imgt input
        hdr.structure <- c("ACCESSION", "GROUP", 
                           "SPECIES", "FUNCTION")
        delim <- "|"
        hdr.structure <- list(header = hdr.structure, delim = delim)    
        id.col <- "GROUP" # use group to identify templates
        if (isolate(input$remove_partial_seqs)) {
            rm.keywords <- c("partial")
        }
    } else {
        # use defined header structure and id column
        # don't use the remove_partial_seqs indicator for other sequences
        hdr.structure <- isolate({header.structure()})
        id.col <- isolate({input$template_header_ID_column})
        rm.duplicated <- isolate(input$remove_duplicated_seqs)
    }
    out <- withWarnings(openPrimeR:::read_templates(seqFile$datapath, 
            hdr.structure = hdr.structure$header, 
            delim = hdr.structure$delim, id.column = id.col, 
            rm.keywords = rm.keywords, remove.duplicates = rm.duplicated,
            gap.character = isolate(gap_char())))
    updateTabsetPanel(session, "main", selected = "Sequences")
    updateTextInput(session, "sample_name", value = seqFile$name) # update analysis identifier
    for (i in seq_along(out$warnings)) {
        warning <- out$warnings[[i]]
        message(warning)
        if (inherits(warning, "TemplateIDColNotFound")) {
            toggleModal(session, "TemplateIDColNotFound")
        }
    }
    for (i in seq_along(out$errors)) {
        error <- out$errors[[i]]
        print(error) # never do warning/message with errors ..
        if (inherits(error, "TemplateHeaderStructure")) {
            toggleModal(session, "TemplateHeaderStructure")
        } else if (inherits(error, "FastaAlphabetError")) {
            toggleModal(session, "FastaAlphabetError")
        } else if (inherits(error, "ID_Column_Not_Found")) {
            toggleModal(session, "IDColumnNotFound")
        } else if (inherits(error, "TemplateFormatIncorrect")) {
            toggleModal(session, "TemplateFormatIncorrect")
        } else {
            toggleModal(session, "TemplateFormatIncorrect")
        }
    }
    if (length(out$errors) != 0) {
        out <- NULL
        session$sendCustomMessage(type = "resetFileInputHandler", 'sequence_file')
    } else {
        out <- out$value
    }
    if (length(out) == 0) {
        rv_templates$SeqTab <- NULL
        validate(need(out, "No sequences available."), errorClass = "critical")
    }
    # switch to template tab and set "all" selector
    updateTabsetPanel(session, "main", selected = "template_view_panel")
    isolate({switch.view.selection("all", input$main, session)})
    return(out)
})
seq.data <- reactive({
    # Loads input template data from seq.data.input(). Resets computed values on input.

    seqs <- seq.data.input() # manual upload of fasta/IMGT DB data
    # update the sliders for uniform allowed regions
    if (length(seqs) != 0) {
        updateSliderInput(session, "uniform_allowed_regions_fw", min = 1, max = max(nchar(seqs$Sequence)))
        updateSliderInput(session, "uniform_allowed_regions_rev", min = 1, max = max(nchar(seqs$Sequence)))
        # activate "confirm" button after template upload
        shinyjs::enable("confirm_uploaded_templates")
    } else {
        shinyjs::disable("confirm_uploaded_templates") # disable the template confirm button
    }
    return(seqs)
})
output$SeqTab <- DT::renderDataTable({
    # Displays the template data in the UI.
    withProgress(message = 'Rendering template table ...', value = 0, {
        data <- switch(input$set_meta_selector,
                "all" = rv_templates$SeqTab, 
                "filtered" = rv_templates$SeqTabFiltered,
                "optimized" = rv_templates$SeqTabOptimized
                )
    })
    validate(need(data, "There is no template data available. Please check your input files and settings."))
    tab <- DT::datatable(asS3(data), caption = "Overview of all uploaded template sequences.", options = list(processing = FALSE), extensions="Responsive")
    return(tab)
})
leader.data.fw <- reactive({
    # fw allowed binding region data

    leaderFile.fw <- rv_cur.input.data$templates_leader
    if (is.null(leaderFile.fw)) {
        return(NULL)
    }
    rm.keywords <- NULL
    if (isolate(rv_templates$load_IMGT_templates)) {
        if (isolate(input$remove_partial_seqs)) {
            rm.keywords <- c("partial")
        }
    }
    gap.char <- isolate(input$gap_char)
    leaders.fw <- withWarnings(openPrimeR:::read.leaders(leaderFile.fw$datapath, "fw", rm.keywords, gap.char))
    for (i in seq_along(leaders.fw$errors)) {
        error <- leaders.fw$errors[[i]]
        print(error)
        if (inherits(error, "FastaAlphabetError")) {
            toggleModal(session, "FastaAlphabetError")
        } else {
            toggleModal(session, "NotifyCouldNotReadFASTA")
        }
    }
    if (length(leaders.fw$errors) != 0) {
        #message("resetting leaders")
        leaders.fw <- NULL
        session$sendCustomMessage(type = "resetFileInputHandler", 'leader_file')
    } else {
        leaders.fw <- leaders.fw$value
    }
    validate(need(length(leaders.fw) != 0, "Could not read the allowed regions for the forward primers."), errorClass = "critical")
    # determine min length before the binding start
    min <- max(leaders.fw$Sequence_Length)
    updateSliderInput(session, "individual_allowed_regions_fw", min = -min, value = c(-min, -1))
    return(leaders.fw)
})
leader.data.rev <- reactive({
    # individual binding region for rev primers
    leaderFile.rev <- input$leader_file_rev
    if (is.null(leaderFile.rev)) {
        return(NULL)
    }
    rm.keywords <- NULL
    if (isolate(rv_templates$load_IMGT_templates)) {
        if (isolate(input$remove_partial_seqs)) {
            rm.keywords <- c("partial")
        }
    }
    gap.char <- isolate(input$gap_char)
    leaders.rev <- withWarnings(openPrimeR:::read.leaders(
        leaderFile.rev$datapath, "rev", rm.keywords, gap.char))
    for (i in seq_along(leaders.rev$errors)) {
        error <- leaders.rev$errors[[i]]
        print(error)
        if (inherits(error, "FastaAlphabetError")) {
            toggleModal(session, "FastaAlphabetError")
        } else {
            toggleModal(session, "NotifyCouldNotReadFASTA")
        }
    }
    if (length(leaders.rev$errors) != 0) {
        session$sendCustomMessage(type = "resetFileInputHandler", 'leader_file_rev')
        leaders.rev <- NULL
    } else {
        leaders.rev <- leaders.rev$value
    }
    validate(need(length(leaders.rev) != 0, "Allowed regions could not be read."), errorClass = "critical")
    max <- 40
    min <- max(leaders.rev$Sequence_Length)
    updateSliderInput(session, "individual_allowed_regions_rev", min = -min, max = max, value = c(-min, -1), step = 1)
    return(leaders.rev)
})

selected.uniform.allowed.regions <- eventReactive(input$uniform_region_confirm_button, { # should be 0 initially
    # when uniform regions have been changed and confirmed, 
    # this updates the reactive values relating to allowed regions. 
    if (input$uniform_region_confirm_button == 0) {
        return(NULL)
    }
    fw <-  input$uniform_allowed_regions_fw
    rev <- input$uniform_allowed_regions_rev
    result <- list("fw" = fw,
               "rev" = rev)
    return(result)
}, ignoreNULL = FALSE) # return NULL in case the button was never pressed

selected.individual.allowed.regions <- eventReactive(input$individual_region_confirm_button, { # should be 0 initially
    # when individual regions have been changed and confirmed, 
    #this updates the reactive values relating to allowed regions. 
    #the button was introduced such that an update of the allowed regions doesn't trigger when dragging the slider.
    if (input$individual_region_confirm_button == 0) {
        return(NULL)
    }
    fw <-  input$individual_allowed_regions_fw
    if (fw[1] == -0.99) { # not adjusted
        fw <- NULL
    }
    rev <- input$individual_allowed_regions_rev
    if (rev[1] == -0.99) { # not adjusted
        rev <- NULL
    }
    result <- list("fw" = fw,
               "rev" = rev)

    return(result)
}, ignoreNULL = FALSE) # trigger initially to return NULL value when confirm button was not pressed.

leader.data.input <- reactive({
# currently loaded allowed binding regions (fw & rev)
    if (length(leader.data.fw()) == 0 && length(leader.data.rev()) == 0 || length(seq.data()) == 0) {
        #message("no leaders loaded")
        return(NULL)
    }
    leaders <- NULL
    leaders <- withWarnings(openPrimeR:::unify.leaders(leader.data.fw(), leader.data.rev(), seq.data(), isolate(gap_char())))
    for (i in seq_along(leaders$errors)) {
        error <- leaders$errors[[i]]
        print(error)
        if (inherits(error, "Leaders_no_matches")) {
            toggleModal(session, "NotifyAllowedNoMatches")
        } else if (!inherits(error, "validation")) {
            toggleModal(session, "UnexpectedError")
        }
    }
    # handle custom warnings with nice pop-ups
    warnings <- leaders$warnings
    for (i in seq_along(warnings)) {
        message(warnings[[i]])
        # don't put up a toggle for not all leaders matched ...
        # this occurs normally when we remove duplicated templates
        #if (inherits(warnings[[i]], "Not_all_leaders_matched")) {
        #    toggleModal(session, "NotifyAllowedNotAllLeadersMatched")
       if (inherits(warnings[[i]], "MissingLeaders")) {
            toggleModal(session, "NotifyAllowedMissing")
       } else if (inherits(warnings[[i]], "RedundantLeaders")) {
            toggleModal(session, "NotifyAllowedRedundant")
        } else if (inherits(warnings[[i]], "LeadersNotFound")) {
            toggleModal(session, "NotifyAllowedNotFound")
        } else if (inherits(warnings[[i]], "AmpliconStartUndefined")) {
            toggleModal(session, "AmpliconStartUndefined")
        }
    }
    if (length(leaders$errors) != 0) {
        message("Setting leaders to NULL because there was an error ...")
        leaders <- NULL
        shinyjs::disable("confirm_uploaded_allowed_regions")
    } else {
        leaders <- leaders$value
        updateTabsetPanel(session, "main", selected = "Sequences")
        # allow pressing the confirm uploaded allowed regions button
        shinyjs::enable("confirm_uploaded_allowed_regions")
    }
    return(leaders)
})
gap_char <- reactive({
    # ensure that gap char is a single character.
    if (length(input$gap_char) == 0 || nchar(input$gap_char) == 0) {
        return("-") # default gap char
    } else {
        return(substring(input$gap_char, 1, 1))
    }
})
leader.data.uniform <- reactive({
    # Loads uniform binding region data using input$uniform_allowed_regions_fw and input$uniform_allowed_regions_rev
    validate(need(seq.data(), "Please specificy the templates first to use the uniform definition of allowed regions."))
    fw.region <- selected.uniform.allowed.regions()$fw
    rev.region <- selected.uniform.allowed.regions()$rev
    if (length(selected.uniform.allowed.regions()) == 0) {
        # customize button wasn't pressed yet -> use the default settings:
        fw.region <- isolate(input$uniform_allowed_regions_fw)
        rev.region <- isolate(input$uniform_allowed_regions_rev)
    }
    leaders <- withWarnings(openPrimeR:::create.uniform.leaders(fw.region, rev.region, seq.data(), isolate(gap_char())))
    for (i in seq_along(leaders$errors)) {
        error <- leaders$errors[[i]]
        print(error)
        toggleModal(session, "UnexpectedError")
    }
    for (i in seq_along(leaders$warnings)) {
        warning <- leaders$warnings[[i]]
        message(warning)
        if (inherits(warning, "AmpliconStartUndefined")) {
            toggleModal(session, "AmpliconStartUndefined")
        }
    }
    if (length(leaders$errors) != 0) {
        leaders <- NULL
    } else {
        leaders <- leaders$value
    }
    return(leaders)
})
leader.data <- reactive({
    # Sets currently active binding region data (either uniform, template-specific, or none)
    if (input$selected_allowed_region_definition == "Uniform") {
        leaders <- leader.data.uniform() # uniform-leaders: template unspecific, generated by positional range info
    } else if (input$selected_allowed_region_definition == "Template-specific") {
        #message("switching to template-specific")
        leaders <- leader.data.input() # template-specific leaders
    } else if (input$selected_allowed_region_definition == "None") {
        leaders <- NULL # no restrictions -> "pure exon" :-)
    }
    return(leaders)
})

LeaderObserver <- observeEvent(c(leader.data(), selected.individual.allowed.regions()), {
    # when leader changes, update primer bindnig region if available
    primer.df <- rv_primers$evaluated_primers
    if (!"primer_coverage" %in% colnames(primer.df)) {
        return()
    }
    template.df <- current.seqs()
    old.template.df <- rv_templates$cvg_all
    # check whether primers correspond to the templates and update coverage
    template.df <- suppressWarnings(try(openPrimeR::update_template_cvg(template.df, primer.df)))
    #print("NEW LEADERS:")
    if (class(template.df) != "try-error") {
        # update the binding positions of the templates relative to the current selected binding region
        primer.df <- openPrimeR:::update_primer_binding_regions(primer.df, template.df, old.template.df)
        # update evaluated_primers with new binding regions
        rv_primers$evaluated_primers <- primer.df
        # update cvg_templates with new binding sites from leader change
        rv_templates$cvg_all <- template.df
    }
})
IMGT_TemplateDataObserver <- observeEvent(input$IMGT_template_button, {
    # retrieves templates from IMGT or local disk (if available) and sets rv_cur.input.data
    # update of partial seqs makes some problems (can only update once)

    # reset current data: necessary to reload data even if filename doesn't change!
    rv_cur.input.data$templates_exon <- NULL
    rv_cur.input.data$templates_leader <- NULL
    rv_cur.input.data$templates_leader_rev <- NULL
    fnames <- retrieve.IMGT.templates(input$IMGT_DB_species, input$IMGT_DB_locus, input$IMGT_DB_function, input$update_IMGT_DB_data, input$remove_partial_seqs)
    if (length(fnames) != 0) {
        rv_templates$load_IMGT_templates <- TRUE # if TRUE -> overwite header structure
        # return fasta filename of exon and leader file in vector
        seqFile <- list("datapath" = fnames[1],  # exon file
                            "name" = basename(fnames[1]))
        leaderFile.fw <- list("datapath" = fnames[2],  # leader file
                            "name" = basename(fnames[2]))
        rv_cur.input.data$templates_exon <- seqFile
        rv_cur.input.data$templates_leader <- leaderFile.fw
        # activate 'confirm templates' button
        shinyjs::enable("IMGT_template_confirm_button")
    } else {
        # disable confirm templates button
        shinyjs::disable("IMGT_template_confirm_button")
        # throw warning:
        session$sendCustomMessage(type='jsCode', list(value = "$('#NotifyIMGT_ConnectionError').modal('show')"))
    }
})

SeqTabObserver <- observe({ 
    # sets the current SeqTab and rv_raw.seqs using the current template data

    #############
    # all data:
    ##############
    if (length(seq.data()) != 0) {
        rv_templates$SeqTab <- view.input.sequences(seq.data())
        rv_templates$raw_seqs <- seq.data()
    }
    if (length(current.seqs()) != 0) { 
        rv_templates$SeqTab <- view.lex.sequences(current.seqs())
        rv_templates$raw_seqs  <- current.seqs()
    }
    if (length(rv_templates$cvg_all) != 0 && length(primer.data()) !=0) {
        rv_templates$SeqTab <- view.cvg.sequences(rv_templates$cvg_all, primer.data())
        rv_templates$raw_seqs <- rv_templates$cvg_all
    }
    #################
    # filtered data:
    ##################
    if (length(rv_templates$cvg_filtered) != 0 && length(current.filtered.primers()) !=0) {
        rv_templates$SeqTabFiltered <- view.cvg.sequences(rv_templates$cvg_filtered, current.filtered.primers())
    }
    ###################
    # optimized data:
    ####################
    if (length(rv_templates$cvg_optimized) != 0 && length(optimal.primers()) != 0) {
        rv_templates$SeqTabOptimized <- view.cvg.sequences(rv_templates$cvg_optimized, optimal.primers())
    }
})
get.exon.data <- reactive({
    # Retrieves the template data with annotations of binding regions.
    if (length(seq.data()) == 0) {
        return(NULL)
    }
    seqs <- seq.data() 
    if (input$selected_allowed_region_definition == "Uniform") { # uniform binding regions selected
        ex.data <- openPrimeR:::add.uniform.leaders.to.seqs(seqs, leader.data()) 
    } else {
        if (length(leader.data()) == 0) { # no individual binding regions available
            #message("no leader data available")
            ex.data <- seqs
        } else { # leader data is available
            #message("loading individual leader data")
            ex.data <- openPrimeR:::get.leader.exon.regions(seqs, leader.data())
            fw.region <- selected.individual.allowed.regions()$fw
            rev.region <- selected.individual.allowed.regions()$rev
            ex.data <- openPrimeR:::adjust_binding_regions(ex.data, fw.region, rev.region)
        }
    }
    validate(need(length(ex.data) != 0, "Could not assign any leader seqs. Check whether your input sequences agree with each other!"), errorClass = "critical")
    # adjust slider for region customization
    if ("Allowed_End_fw_initial" %in% colnames(ex.data)) {
        # adjust max of fw after template.df has been determined ...
        max <- max(nchar(ex.data$Sequence) - ex.data$Allowed_End_fw_initial) - 1
        updateSliderInput(session, "individual_allowed_regions_fw", max = max)
    }
    if ("Allowed_Start_rev_initial" %in% colnames(ex.data)) {
        # adjust max of rev allowed region
        max <- max(ex.data$Allowed_Start_rev_initial) - 1
        updateSliderInput(session, "individual_allowed_regions_rev", max = max)
    }
    return(ex.data)
})

current.seqs <- reactive({
    # The currently loaded templates, with target binding regions if available.

    lex.seqs <- get.exon.data() # seqs with binding annotations
    seqs <- seq.data() # seqs without binding region annotations
    if (length(lex.seqs) != 0) {
        seqs <- lex.seqs
    }
    # update binding region by secondary structure
    opti.regions <- optimized.regions.structure()$Intervals
    if (length(seqs) != 0 && length(opti.regions) != 0 && length(opti.regions[[1]]) == nrow(seqs)) {
        # update optimized regions only for the matching template set
        seqs <- openPrimeR:::update.binding.regions(seqs, opti.regions)
    }
    new.seqs <- optimized.regions.conservation() # doesn't use the current seqs ...
    # update binding region by conservation
    #new.seqs <- openPrimeR:::select_regions_by_conservation(seqs
                                    #gap.char = gap_char(),
                                    #win.len = 40, by.group = TRUE,
                                    #direction = input$design_drection)
    if (length(new.seqs) != 0) {
        seqs <- new.seqs
    }
    return(seqs)
})
header.structure <- reactive({
# Determines the fasta header structure info required for loading the templates.
    delim <- input$template_header_delim
    hdr.str <- input$template_header_structure
    # in case settings have not been shown yet, use defaults:
    if (length(delim) == 0) {
        delim <- "|"
    }
    if (length(hdr.str) == 0) {
        # hdr.structure not defined yet in the UI
        # only use accession to be able to load any fasta file
        hdr.str <- c("ACCESSION") 
    }
    if (length(hdr.str) == 0 && length(delim) == 0) {
        delim <- ""
        hdr.str <- ""
    } else if (length(hdr.str) == 0) {
        hdr.str <- ""            
    }
    return(list("header" = hdr.str, "delim" = delim))
})
conservation_plot_height <- reactive({
    if (length(current.seqs()) == 0) {
        return(800)
    }
    nbr.groups <- length(unique(current.seqs()$Group))
    height <- openPrimeR:::get.plot.height(ceiling(nbr.groups/2), 300, 600)
})

output$header_structure <- renderUI({
# Displays the available template header fields to the user in the frontend for selection.
  fields <- c("ACCESSION", "GROUP", "SPECIES", "FUNCTION")
  ## using selectizeInput with drag_drop and DT
  selectizeInput("template_header_structure", tagList(icon("list-alt", lib = "glyphicon"), "FASTA header fields"), choices  = fields,
   selected = "ACCESSION", multiple = TRUE,
   options = list(plugins = list('remove_button', 'drag_drop')))
})
output$template_conservation_plot <- renderPlot({
    validate(need(attr(current.seqs(), "entropies"), "Conservation is not available yet."))
    entropies <- attr(current.seqs(), "entropies")
    alignments <- attr(current.seqs(), "alignments")
    openPrimeR:::plot_conservation(entropies, alignments, current.seqs(), 
                                   gap.char = isolate(gap_char()))
}, height = conservation_plot_height)


output$template_secondary_plot <- renderPlot({
    validate(need(optimized.regions.structure(), "Template binding regions haven't been optimized yet."))
    fold.df <- optimized.regions.structure()$Foldings
    openPrimeR:::plot_template_structure(fold.df)
})
optimized.regions.structure <- eventReactive(input$modify_binding_regions_secondary_structures, {
    # modify template target region based on secondary structure
    if (input$modify_binding_regions_secondary_structures == 0 || length(get.exon.data()) == 0) { # no target region annotation available ...
        return(NULL)
    }
    isolate(annealing.temp <- annealing.temperature()) # don't trigger
    result <- openPrimeR:::optimize.template.binding.regions.dir(get.exon.data(), annealing.temp, input$allowed_primer_length, input$design_direction)
    # result: consists of 'Intervals' (new binding regions) and 'Foldings' (data frame with DeltaG information)
    if (length(result) == 0) { # nothing could be changed (no regions defined)
        return()
    }
    updateTabsetPanel(session, "main", selected = "template_view_panel") # update view to template tab
    return(result)
}, ignoreNULL = FALSE) # trigger also on NULL to return something

optimized.regions.conservation <- eventReactive(input$modify_binding_regions_conservation, {
    # modify template target region based on secondary structure
    if (input$modify_binding_regions_conservation == 0 || length(get.exon.data()) == 0) { # no target region annotation available ...
        return(NULL)
    }
    result <- openPrimeR:::select_regions_by_conservation(get.exon.data(), 
                                    gap.char = isolate(gap_char()),
                                    win.len = 40, by.group = TRUE,
                                    direction = input$design_drection)
    if (length(result) == 0) { # nothing could be changed (no regions defined)
        return()
    }
    # update view to template tab
    updateTabsetPanel(session, "main", selected = "template_view_panel") 
    return(result)
}, ignoreNULL = FALSE) # trigger also on NULL to return something

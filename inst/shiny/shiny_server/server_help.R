######## Shiny server functionalities for the help pages

#HelpObserver <- observe({
    ## HelpObserver is just there to summarize all help events input help
    ##############
    # TEMPLATES
    ###############
    observeEvent(input$help_input_templates, {
        view.template.help.element(session, "help_input_templates_overview")
    })
    observeEvent(input$help_input_templates_comparison, {
        view.template.help.element(session, "help_input_templates_overview")
    })
    observeEvent(input$help_input_templates_uniform, {
        view.template.help.element(session, "help_input_templates_allowed")
    })
    observeEvent(input$help_input_templates_allowed, {
        view.template.help.element(session, "help_input_templates_allowed")
    })
    observeEvent(input$help_input_templates_header, {
        view.template.help.element(session, "help_input_templates_header")
    })
    ################
    # PRIMERS
    ##################
    observeEvent(input$help_input_primers_comparison, {
        # Comparison primer input
        view.primer.help.element(session, "help_input_primers_comparison")
    })
    observeEvent(input$help_input_primers_overview, {
        # Input of primers for analysis
        view.primer.help.element(session, "help_input_primers_overview")
    })
    observeEvent(input$help_init_initialization, {
        # specifying the properties of design primers
        view.opti.help.element(session, "help_init_initialization")
    })
    #############
    # SETTINGS
    ##############
    ########
    # a) OVERALL SETTINGS
    #########
    observeEvent(input$help_settings_overview, {
        # overall settings info
        view.constraint.general.help("help_constraints_overview", session)
    })
    observeEvent(input$help_coverage_conditions, {
        # coverage conditions
        view.cvg.help("help_tab_coverage_basic", session)
    })
    observeEvent(input$help_pcr_settings, {
        # PCR conditions
        view.PCR.help(session)
    })
    #######
    # b) CONSTRAINTS
    #######
    # help for constraints
    observeEvent(input$help_overview_filters, {
        view.filter.help("help_tab_overview_filters", session)
    })
    observeEvent(input$help_primer_coverage, {
        view.filter.help("help_tab_primer_coverage", session)
    })
    observeEvent(input$help_primer_length, {
        view.filter.help("help_tab_primer_length", session)
    })
    observeEvent(input$help_gc_clamp, {
        view.filter.help("help_tab_gc_clamp", session)
    })
    observeEvent(input$help_gc_ratio, {
        view.filter.help("help_tab_gc_ratio", session)
    })
    observeEvent(input$help_run_length, {
        view.filter.help("help_tab_run_length", session)
    })
    observeEvent(input$help_repeat_length, {
        view.filter.help("help_tab_repeat_length", session)
    })
    observeEvent(input$help_melting_temperature, {
        view.filter.help("help_tab_melting_temperature", session)
    })
    observeEvent(input$help_opti_melting_temp, {
        # melting temp diff
        view.filter.help("help_tab_melting_temperature_diff", session)
    })
    observeEvent(input$help_secondary_structure, {
        view.filter.help("help_tab_secondary_structure", session)
    })
    observeEvent(input$help_primer_specificity, {
        view.filter.help("help_tab_primer_specificity", session)
    })
    observeEvent(input$help_cross_complementarity, {
        view.filter.help("help_tab_cross_complementarity", session)
    })
    observeEvent(input$help_self_complementarity, {
        view.filter.help("help_tab_self_complementarity", session)
    })
    ########
    # c) Coverage conditions
    #########
    observeEvent(input$help_primer_efficiency, {
        view.cvg.help("help_tab_primer_efficiency", session)
    })
    observeEvent(input$help_coverage_model, {
        view.cvg.help("help_tab_coverage_model", session)
    })
    observeEvent(input$help_annealing_DeltaG, {
        view.cvg.help("help_tab_annealing_DeltaG", session)
    })
    observeEvent(input$help_codon_design, {
        view.cvg.help("help_tab_codon_design", session)
    })
    ################
    # ANALYSIS
    #################
    # EVALUATION
    observeEvent(input$help_eval_eval, {
        view.eval.help.entry(session, "help_eval_eval")
    })
    # COMPARISON
    observeEvent(input$help_compare, {
        view.compare.help(session)
    })
    # OPTIMIZATION
    observeEvent(input$help_opti_optimization, {
        view.opti.help.element(session, "opti_optimization_help")
    })
    observeEvent(input$help_opti_optimization_overview, {
        view.opti.help.element(session, "opti_optimization_help")
    })
    observeEvent(input$help_modify_binding_regions_secondary_structures, {
        # optimize template secondary structures
        view.opti.help.element(session, "opti_templates_help")
    })
######### help pages end
#})
############################ 

output$helpPage <- renderUI({
    # read help page from external file and output it as HTML
    help.file <- "help.html"
    content <- readChar(help.file, file.info(help.file)$size)
    return(HTML(content))
})

output$iupac_ambig_table <- renderTable({
    # output a table with IUPAC ambiguities conversions
    iupac.table <- data.frame("IUPAC Code" = names(Biostrings::IUPAC_CODE_MAP), "Nucleotide" = sapply(strsplit(Biostrings::IUPAC_CODE_MAP, 
    split = ""), function(x) paste(x, collapse = " or ")))
    colnames(iupac.table) = c("IUPAC Code", "Nucleotide")
    return(iupac.table)
}, include.rownames = FALSE, caption = "Table 1: IUPAC ambiguity codes", caption.placement = getOption("xtable.caption.placement", 
    "bottom"), caption.width = getOption("xtable.caption.width", NULL))

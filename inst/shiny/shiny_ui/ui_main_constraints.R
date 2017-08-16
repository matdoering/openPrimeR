########
# Constraints main panel in UI
#######
tabPanel("Constraints",
    value = "constraints_view",
    icon  = icon("resize-small", lib = "glyphicon"),
    br(),
    selectizeInput("selected_constraint_result",  # select a constraint result to display
            "Selected result",
             choices = list(
                        Primers = 
                            c(`Constraint Overview` = "overview",
                              `Dimerization` = "dimerization",
                              `Other` = "other_constraints"),
                        Coverage = 
                            c(`Mismatches` = "no_mismatches",
                              `Other` = "other_cvg_constraints"),
                              #`Primer Efficiencies` = "primer_efficiency"),
                        Templates = 
                            c(`Secondary Structures` = "template_secondary",
                              `Sequence Conservation` = "template_conservation")
            )
    ),
    #######
    # Overview of constraints:
    ###########
    conditionalPanel("input.set_meta_selector != 'all' &&
                        input.selected_constraint_result == 'overview'", 
            # show filtering/optimization stats
            selectInput("selected_filtering_plot", "Output", 
                               c("Evaluation" = "evaluation",
                                 "Filtering overview" = "overview", 
                                 "Filtering coverage" = "coverage", 
                                 "Filtering runtime" = "runtime")
            ),
            conditionalPanel(
                # filtering overview
                condition = "input.selected_filtering_plot == 'overview'", 
                fluidRow(
                    column(12, plotOutput("filtering_stats"), align= "center")
               )
             ),
             conditionalPanel(
                # filtering coverage
                condition = "input.selected_filtering_plot == 'coverage'", 
                fluidRow(
                    column(12, plotOutput("filtering_stats_cvg"), align = "center"),
                    column(12, plotOutput("exclusion_stats"), align = "center")
                )
             ),
             conditionalPanel(
                # filtering runtime
                condition = "input.selected_filtering_plot == 'runtime'", 
                fluidRow(
                    column(12, plotOutput("filtering_runtime"), align="center")
                )
            )
        ),
    conditionalPanel("input.selected_constraint_result == 'overview'",
        conditionalPanel("input.set_meta_selector == 'all' ||
                          input.selected_filtering_plot == 'evaluation'",  
            # evaluation stats
            #column(12, plotOutput("constraint_fulfillment_plot"), align="center"), # removed due to redundancy with 'constraint_stats' plot which provides more infos
            # text summary of constraints
            column(12, uiOutput("ConstraintsTotal")),
            br(),
            # fulfillment matrix:
            column(12, plotOutput("constraint_stats"), align="center"),
            # deviation of constraints:
            column(12, plotOutput("constraint_deviations"), align = "center")
        )
    ),
    ########
    # Dimerization
    ##########
        conditionalPanel("input.selected_constraint_result == 'dimerization'",
            # dimerization results
            selectizeInput("selected_dimerization_data", "Dimerization type", c("Self-Dimerization", "Cross-Dimerization"),
            options = list(
                placeholder = 'Please select a type of dimerization.',
                onInitialize = I('function() { this.setValue(""); }')
                )
            ), # selectize: don't select anything without user input
            selectInput("selected_dimerization_result", "Selected result", c("Summary" = "summary", 
                                "Details" = "details"
            )),
            conditionalPanel(
                # dimerization summary
                condition = "input.selected_dimerization_result == 'summary'", 
                htmlOutput("dimer_text"),
                plotOutput("dimer_distribution"),
                DT::dataTableOutput("dimer_table") # warning (structure()) due to NULL value (caused by shiny, WONTFIX)
             ),
             conditionalPanel(
                # dimerization details
                condition = "input.selected_dimerization_result == 'details'", 
                    DT::dataTableOutput("dimer_data")
                    #bsModal("DimerDetail", "All possible dimerizations", "", size = "large",
                        #DT::dataTableOutput("dimer_detail"))
            )
        ),
    ########
    # Other constraints
    ###########
    conditionalPanel("input.selected_constraint_result == 'other_constraints'",
        selectizeInput("selected_other_result", 
            "Selected constraints",
            multiple = TRUE,
            choices = NULL, # set choices dynamically based on the settings
             options = list(
                placeholder = 'Please select some constraints',
                onInitialize = I('function() { this.setValue(""); }')
             )
            #choices = c("Primer coverage" = "primer_coverage",
                        #"Primer length" = "primer_length",
                        #"Melting temperature range" = "melting_temp_range",
                        #"Melting temperature differences" = "melting_temp_diff",
                        #"GC clamp" = "gc_clamp",
                        #"Number of runs" = "no_runs",
                        #"Number of repeats" = "no_repeats",
                        #"Secondary structures" = "secondary_structure",
                        #"Primer specificity" = "primer_specificity",
                        #"Self-Dimerization" = "self_dimerization",
                        #"Cross-Dimerization" = "cross_dimerization")
        ),
        plotOutput("constraint_plot_histogram")
    ),

    ##########
    # Coverage constraints
    ##########
    conditionalPanel("input.selected_constraint_result == 'no_mismatches'",
        # nbr mismatch plot
        plotOutput("constraint_plots_no_mismatches"),
        selectInput("selected_primer_set_mismatches_direction", "Selected primer direction", c("fw", "rev")),
        br(),
        DT::dataTableOutput('mismatch_table')
    ),
    # other cvg constraints:
    conditionalPanel("input.selected_constraint_result == 'other_cvg_constraints'",
        selectizeInput("selected_cvg_constraints", 
            "Selected coverage constraints",
            multiple = TRUE,
            choices = NULL, # set choices dynamically based on the selected settings
            options = list(
                placeholder = 'Please select a coverage constraint',
                onInitialize = I('function() { this.setValue(""); }')
                )
            #choices = c("Primer efficiency" = "primer_efficiency",
                        #"Primer annealing" = "annealing_DeltaG",
                        #"Stop codons" = "stop_codon",
                        #"Terminal mismatches" = "terminal_mismatch_pos",
                        #"Coverage Model FPR" = "coverage_model")
        ),
        plotOutput("constraint_plots_cvg_constraints")
    ),
    ########
    # TEMPLATE CONSTRAINTS
    ##########
    conditionalPanel("input.selected_constraint_result == 'template_secondary'",
        plotOutput("template_secondary_plot")
    ),
    conditionalPanel("input.selected_constraint_result == 'template_conservation'",
        plotOutput("template_conservation_plot")
    )

)


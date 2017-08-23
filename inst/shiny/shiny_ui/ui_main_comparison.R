######
# Comparison main panel in UI
######
tabPanel("Comparison",
    icon = icon("scale", lib = "glyphicon"),
    br(),
    ###################
    selectizeInput("selected_comparison_plot",  # comparison result to display
                "Selected result", 
              choices = list(
              Data = c(`Loaded Data` = "loaded_data"),
              Coverage = c(`Coverage` = "coverage_overview", 
                `Binding Regions` = "binding_regions",
                `Coverage vs Size` = "coverage_vs_size",
                `Mismatches` = "mismatches",
                `Coverage Constraints` = "cvg_constraints"),
              Constraints = c(`Constraint Fulfillment` = "constraints_overview",
                `Constraint Table`  = "constraints_table",
                `Constraint Details` = "other_constraints")
    )),
    conditionalPanel(
        condition = "input.selected_comparison_plot == 'loaded_data'",  # overview of loaded data
        fluidRow(
            column(12, 
                DT::dataTableOutput("uploaded_comp_data")
             )
        )
    ),
    conditionalPanel(
        condition = "input.selected_comparison_plot == 'constraints_overview'",  # overview of constraints
            fluidRow(
                column(12, 
                    p("Modifying the selected constraints will update the plot."),
                    uiOutput("comparison_plot_evaluation_ui")
                ),
                column(12,
                    plotOutput("comparison_plot_deviation")
                )
            )
     ),
     conditionalPanel(
        condition = "input.selected_comparison_plot == 'constraints_table'",  # comparison of primer set constraint stats
        fluidRow(
            column(12, 
                DT::dataTableOutput('comparison_overview_table'),
                align= "center"
            )
        )
     ),
     conditionalPanel(
        condition = "input.selected_comparison_plot == 'coverage_overview'",  # cvg comparison
        fluidRow(
            column(12, 
                plotOutput("comparison_plot_cvg"),align="center"
            ),
            column(12, 
                DT::dataTableOutput("comparison_stats"), align= "center"
            )
        )
     ),
     conditionalPanel(
        condition = "input.selected_comparison_plot == 'coverage_vs_size'",  # cvg vs set size
        fluidRow(
            column(12, 
                plotOutput("cvg_vs_size_plot"),align="center"
            )
        )
     ),
     conditionalPanel(
        # coverage constraints
        condition = "input.selected_comparison_plot == 'cvg_constraints'",
        selectizeInput("selected_cvg_comp_constraints", 
            "Selected coverage constraints",
            multiple = TRUE,
            choices = NULL,
            options = list(
                placeholder = 'Please select a coverage constraint',
                onInitialize = I('function() { this.setValue(""); }')
             )
            #choices = c("Primer efficiency" = "primer_efficiency",
                        #"Primer annealing" = "annealing_DeltaG",
                        #"Coverage model FPR" = "coverage_model",
                        #"Stop codons" = "stop_codon",
                        #"Terminal mismatches" = "terminal_mismatch_pos")
        ),
        fluidRow(
            column(12, 
                plotOutput("comp_cvg_constraints"),align="center"
            )
        )
     ),

     conditionalPanel(
        condition = "input.selected_comparison_plot == 'other_constraints'", 
         selectizeInput("selected_other_plot",  # comparison result to display
                "Selected constraints", 
                 choices = NULL,
                 options = list(
                    placeholder = 'Please select some constraints',
                    onInitialize = I('function() { this.setValue(""); }')
                 ),
                multiple = TRUE
                #c(
                #"Primer coverage" = "primer_coverage",
                #"Self-Dimerization" = "self_dimerization",
                #"Cross-Dimerization" = "cross_dimerization",
                #"Primer length" = "primer_length",
                #"Melting temperature range" = "melting_temp_range",
                #"Melting temperature differences" = "melting_temp_diff",
                #"GC clamp" = "gc_clamp",
                #"GC ratio" = "gc_ratio",
                #"Number of runs" = "no_runs",
                #"Number of repeats" = "no_repeats",
                #"Secondary structures" = "secondary_structure",
                #"Primer specificity" = "primer_specificity"),
        ),
        plotOutput("comparison_plot_constraint")
    ),
    conditionalPanel(
        condition = "input.selected_comparison_plot == 'binding_regions'", # binding region comparison
        selectInput("primer_comparison_relation", "Plot relative to which binding region?", c("Forward" = "fw", "Reverse" = "rev")),
        fluidRow(
            column(12, 
                plotOutput("comparison_plot_regions"),align= "center"
            )
        )
     ),
    conditionalPanel(
        "input.selected_comparison_plot == 'mismatches'",  # mismatch comparison
        fluidRow(
            column(12, 
                plotOutput("comparison_plot_mismatches"),align= "center"

            )
        )
    )
)


######
# Coverage main panel in UI
######
tabPanel("Coverage", 
    icon = icon("stats", lib = "glyphicon"),
    value = "coverage_view_panel",
    br(),
    fluidRow(
        column(4,
        selectizeInput("selected_cvg_plot",  # select specific result
            "Selected result", 
            choices = list(
            Templates = 
                c(`Template Coverage` = "overview", 
                  `Mismatch Template Coverage` = "overview_mismatches"),
            Primers = 
                c(`Primer Coverage` = "individual_coverage", 
                `Mismatch Primer Coverage` = "individual_coverage_mismatches",
                `Subset Coverage` = "subset_coverage"),
            Binding = 
                c(`Binding Regions` = "primer_binding",
                `Primer View` = "primer_view")
        ))),
        column(4,
            selectizeInput( # select gene group
                "selected_group_coverage", "Selected group",
                choices = c("all"), multiple = TRUE, selected = "all"
        ))
    ),
    conditionalPanel(
        # coverage overview for templates
        condition = "input.selected_cvg_plot == 'overview'", 
        fluidRow(
            htmlOutput("CoverageTotal"),
            br(),
            column(3, wellPanel(uiOutput("CoverageBox"))),
            column(9, uiOutput("coverage_primer_per_group_ui"), align = "center"),
            selectizeInput( # select coverage definition for stat table
                "selected_cvg_def_stats", "Selected coverage definition",
                choices = c("Expected Coverage" = "Expected_Coverage", "Identity Coverage" = "Identity_Coverage"), 
                multiple = FALSE, selected = "expected_cvg"
            ),
            column(12, DT::dataTableOutput("cvg_stats"), align = "center")
        )
     ),
     conditionalPanel(
        # mismatch coverage overview for templates
        condition = "input.selected_cvg_plot == 'overview_mismatches'", 
        fluidRow(
            column(12, uiOutput("template_coverage_mismatch_ui"), align = "center")
            #column(12, plotOutput("template_coverage_mismatch"), align = "center")

        ),
        fluidRow(
            column(12,  
            selectizeInput( # select nbr of allowed mismatches for table
                "allowed_mm_cvg_stats", "Allowed number of mismatches",
                # TODO
                NULL, 
                options = list(
                    placeholder = 'Select a number of mismatches',
                    onInitialize = I('function() { this.setValue(""); }')
                )
            )),
             span(paste("Retrieve template coverage statistics",
                "for primers binding with the allowed number of mismatches.")
            ),
            br(),
            column(12, DT::dataTableOutput("cvg_stats_mismatches"), align = "center"),
            br(), br(), br(), br(), br(),br(),br(),br(), br() # extra space for showing options of mismatch number selector
       )
     ),
     conditionalPanel(
        # individual primer coverage (for defined cvg condition)
        condition = "input.selected_cvg_plot == 'individual_coverage'", 
        fluidRow(
            column(12, plotOutput("Coverage_Primer"), align = "center")
        )
     ),
     conditionalPanel(
        # individual primer coverage (for all mismatch settings)
        condition = "input.selected_cvg_plot == 'individual_coverage_mismatches'", 
        fluidRow(
            # coverage stats wrt mismatches:
            column(12, DT::dataTableOutput("primer_cvg_stats_mismatch"), align = "center"),
            # primer cvg wrt mismatches:
            column(12, plotOutput("Coverage_Primer_mismatches"), align = "center")
        )
    ),
     conditionalPanel(
        # primer binding plot
        condition = "input.selected_cvg_plot == 'primer_binding'", 
        selectInput("primer_location_plot_rel", "Relation to binding region", c("fw", "rev")),
        fluidRow(
            column(10, plotOutput("primer_binding_regions"), align = "center")
        )
     ),
     conditionalPanel(
        # individual primer binding sites
        condition = "input.selected_cvg_plot == 'primer_view'", 
            fluidRow(
            column(4, 
                selectInput("primer_plot_rel", "Relation to binding region", c("fw", "rev"))),
            column(4,
                selectInput("selected_primer", "Selected primer","", multiple = TRUE))
            ),
            plotOutput("primer_plot")
     ),
     conditionalPanel(
        # coverage of primer subsets and subset selection
        condition = "input.selected_cvg_plot == 'subset_coverage'", 
        fluidRow(
            column(12, 
                plotOutput("primer_subset_coverage"),
                align = "center"
            )
        ),
        fluidRow(
            column(12,
                selectizeInput("selected_subset_size", "Selected subset", NULL,
                options = list(
                    placeholder = 'Please select a size of the primer subset.',
                    onInitialize = I('function() { this.setValue(""); }')
                )
            ),
            span(paste("Retrieve a subset of primers maximizing",
                "the coverage without checking for any other constraints.")
            ),
            br(),
            DT::dataTableOutput('primer_subset_table'), 
            br(), br(), br(), br(), br(),br(),br(),br(), br() # extra space for size selector
        )
    ))
)


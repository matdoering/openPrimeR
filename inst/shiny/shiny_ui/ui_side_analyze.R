########
# UI side panel for analysis of primers
########
tabPanel("Analyze", 
    icon = icon("cogs"),
    value = "analyze_panel",
    br(),
    h2(class="inline", "Analyze"),
    br(),
    br(),
    conditionalPanel("input.primer_analysis_type == 'evaluate'",
        ######
        # EVALUATION PANEL
        ######
        div(p("Evaluate the physicochemical properties of the input primers.",
             create.help.button("eval_eval")),
             class = "one"
        ),
        actionButton("evaluateButton", 
            "Compute primer properties", 
            icon = icon("bar-chart"), class="actionStyleRun btn-primary"
        ),
        bsTooltip("evaluateButton", 
            "Compute the currently active primer properties but do not perform any filtering.",
            "right", options = list(container = "body")
        ),
        br(),br(),
        actionButton("quickFilterButton", 
            "Filter primers", 
            icon = icon("filter", lib = "glyphicon"), class="actionStyleRun btn-primary"
        ),
        bsTooltip("quickFilterButton", 
            "Determine the set of primers fulfilling the active primer selection criteria.",
            "right", options = list(container = "body")
        ),
        br(),br(),
        # adapter check button
        actionButton("check_adapters", 
                     "Check restriction sites",
                      icon = icon("scissors"), class="actionStyleRun btn-primary"
        ),
        bsTooltip("check_adapters", 
            "Identify whether the input primers might contain restriction sites.",
                "right", options = list(container = "body"))
        ), 
        conditionalPanel("input.primer_analysis_type == 'design'",
        ##########
        # DESIGN PANEL
        ##########
        div(p("Design a new set of primers for the loaded template sequences using the currently selected constraints as a selection criterion.", 
            create.help.button("opti_optimization_overview")),
            class = "one"
        ),
        #bsCollapse(id = "design_options_algorithms_collapse", 
        bsCollapse(id = "design_options_algorithms_collapse", 
            open = "design_algo_panel",
            bsCollapsePanel(tagList(
                 icon("menu-hamburger", lib = "glyphicon"), 
                 "Design options"),
                 value = "design_algo_panel",
                div(p("The optimization algorithm for primer design influences two parameters: the size of the constructed primer set and the runtime. Integer linear programs (ILPs) gurantee minimal primer sets, but their worst-case runtime is exponential. Greedy algorithms, are faster (O n log n), but may yield slightly larger primer sets. ILPs are recommended when the size of the designed primer set should be as small as feasible or when primers are to be designed for few templates, or when the runtime is of no major concern.", create.help.button("opti_algo_ref")), class = "two"),
                 bsTooltip("help_opti_algo_ref", 
                        HTML("References:<br>Pearson et al. (1996)<br>Bashir et al. (2007)"),
                        "right", options = list(container = "body"), 
                        trigger = "click"
                ),
                selectInput("optimization_algorithm", tagList(tagList(icon("cogs"), "Optimization algorithm")), choices  = c("Greedy", "ILP"), selected = "Greedy"), # used optimization algorithm
                bsTooltip("optimization_algorithm", 
                    "The algorithm for maximizing the coverage while minimizing the size of the primer set.",
                        "right", options = list(container = "body"), 
                        trigger = "click"
                ),
                #############
                # expert design options:
                #############
                bsCollapse(id = "design_expert_options",
                open = NULL,
                bsCollapsePanel(
                    tagList(
                        icon("menu-hamburger", lib = "glyphicon"), 
                        "Expert options"
                    ),
                    value = "design_expert_panel",
                checkboxInput("relaxation_active", 
                    "Constraint relaxation",
                    TRUE
                ),
                bsTooltip("relaxation_active", 
                        HTML("If the target sequence coverage cannot be reached, <br>the constraints are adjusted until the target coverage can be reached."),
                        "right", options = list(container = "body")
                ),
                conditionalPanel("input.relaxation_active", 
                    sliderInput("required_opti_cvg", 
                        tagList("Target coverage"), 
                        min = 0, max = 1, value = c(min = 1)
                    ), 
                    bsTooltip("required_opti_cvg", 
                        HTML("The percentage of template sequences <br>that should be covered by the designed primers."),
                        "right", options = list(container = "body")
                    )
                ) # conditional panel for relaxation closes
                )), # design expert panel closes
                # design options next button:
                div(class = "rightAligned", 
                    actionButton("confirm_algorithm_choice", 
                        "Confirm design options", icon = icon("check"), 
                        class = "actionStyleSmall",
                        style="color: #fff; background-color: #337ab7; border-color: #2e6da4"))
        ), # close design options panel
        #########
        # TEMPLATE REGION PANEL TODO
        #########
        bsCollapsePanel(tagList(
                 icon("book", lib = "glyphicon"), 
                 "Optimization of binding regions"),
                 value = "design_optimize_binding",
                 div(p("The selected primer binding regions in the template sequences may not be well-suited for designing primers in case that they are subject to the formation of secondary structures of if they are not highly conserved. In these cases, you can obtain better results if you narrow down the possible primer binding range by considering these two factors."),
            class = "two"
                ),
        # optimization of binding region: template secondary structure 
            actionButton("modify_binding_regions_secondary_structures", "Avoid secondary structures", icon = icon("eraser"), class = "actionStyleRun btn-primary"), 
            bsTooltip("modify_binding_regions_secondary_structures", 
                 "Adjust the primer binding regions such that regions exhibiting secondary structures  are avoided.",
                "right", options = list(container = "body")
            ),
            br(),br(),
            actionButton("modify_binding_regions_conservation", "Select conserved binding regions", icon = icon("signal"), class = "actionStyleRun btn-primary"), 
            bsTooltip("modify_binding_regions_conservation", 
                 "Adjust the primer binding regions such that the most conserved regions are selected.",
                "right", options = list(container = "body")
            )
        ) # close panel
        ), # close collapse
        ##########
        # design button
        #######
        actionButton("designButton", "Design primers", 
                icon = icon("sort-amount-desc"),
                class="actionStyleRun btn-primary"),
        bsTooltip("designButton", 
            "Find a minimal set of primers maximizing the coverage subject to the specified design constraints.",
            "right", options = list(container = "body")
        )
    ), 
    conditionalPanel("input.primer_analysis_type == 'compare'",
        ############
        # COMPARISON PANEL
        ############
        div(p("Compare the properties of several primer sets with each other.",
            create.help.button("compare")),
            class = "one"
        ),
        # comparison button
        actionButton("compare_primers", "Compare primer sets", icon = icon("scale", lib = "glyphicon"), class = "actionStyleRun btn-primary"),
        bsTooltip("compare_primers", 
            "Compare the properties of the selected primer sets with each other.",
            "right", options = list(container = "body")
        ),
        br(),
        br(),
        # filter button for comparison
        actionButton("quickFilterButton_compare", "Filter primer sets", icon = icon("filter", lib = "glyphicon"), class="actionStyleRun btn-primary"),br(),
        bsTooltip("quickFilterButton_compare", 
            "Filter all selected primer sets according to the currently active constraints.",
            "right", options = list(container = "body")
        )
    ) 
)

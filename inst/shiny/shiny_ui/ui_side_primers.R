####
# UI side panel for primer input
######
tabPanel("Primers",
    value = "primer_input_tab",
    icon = icon("tag", lib = "glyphicon"),
    br(),                    
    h2(class="inline", "Primers"),
    br(),
    br(),
    conditionalPanel("input.primer_analysis_type == 'evaluate'",
        # load input primers if we're not designing new primers
        withTags(
            div(class = "one", 
            p("For immunological applications, we provide a number 
               of published primer sets that can be loaded directly. 
               For other applications please input a FASTA or CSV file 
               containing the primers to be analyzed.", 
            create.help.button("input_primers_overview"))
        )),
        bsTooltip("help_input_primers", 
            "View help on the input of primers.",
            "right", options = list(container = "body")
        )
    ),
    conditionalPanel("input.primer_analysis_type == 'evaluate'", 
        ###########
        # INPUT FOR PRIMER EVALUATION
        ###########
        # treatment of ambiguities
        radioButtons("use_ambig", 
            "Treatment of IUPAC ambiguity codes", 
            c("None" = "none", "Merge" = "merge", "Disambiguate" = "unmerge"), 
            selected = "none", inline = TRUE),
        bsTooltip("use_ambig", 
            paste("Whether similar primers in the input shall be left as they are",
            "(\\'None\\'), merged using ambiguity codes (\\'Merge\\'),",
            "or disambiguated (\\'Disambiguate\\')."),
        "right", options = list(container = "body")),
        ##################
        # Selection of primer source
        #####################

        # load personal primers or supplied primers?
        radioButtons("primer_upload_choice",
            tagList(icon("floppy-disk", lib = "glyphicon"), "Primer source"),
            choices = c("Supplied primers" = "available_primers",
                        "Personal primers" = "personal_primers"), 
            selected = "available_primers", inline=TRUE
        ),
        bsTooltip("primer_upload_choice", 
            "Whether you want to load supplied primers or your own primers.",
            "right", options = list(container = "body")
        ),

        # checkbox: load evaluated primer csv or fasta? 
        checkboxInput("load_eval_primers", 
                    "Load evaluated primers (CSV)", 
                    TRUE),
        bsTooltip("load_eval_primers", 
                    "Whether pre-evaluated primers (CSV) or raw primers (FASTA) shall be loaded.",
                    "right", options = list(container = "body")
        ),
        conditionalPanel("input.primer_upload_choice == 'available_primers'",
            # Load Supplied IMGT primers
            selectizeInput("IMGT_primers", 
                tagList(icon = icon("tag", lib = "glyphicon"),
                    "Available primers"),
                choices = NULL, selected = NULL,
                options = list(
                    placeholder = 'Please select one of the available primer data sets',
                    onInitialize = I('function() { this.setValue(""); }')
                )
            )
        ),
        conditionalPanel("input.primer_upload_choice == 'personal_primers'", 
            ########
            # Load personal primers
            ######
            bsCollapse(id = "primer_upload_collapse",
                # primer upload panel
                open = "primer_header_structure_panel",
                    #######
                    # INPUT OPTIONS FOR PRIMERS
                    #############
                    bsCollapsePanel(
                    # header structure panel: keywords for fw/rev primers
                    tagList(icon("menu-hamburger", lib = "glyphicon"), 
                    "Primer input options"),
                    value = "primer_header_structure_panel",
                    style = "primary",
                    ##########
                    # only show FASTA options if loading non-evaluated primers
                    ########
                    conditionalPanel("!input.load_eval_primers", 
                        div(p("Please specify the keywords in the headers 
                           of the FASTA file that identify the primer 
                           directionalities."), class="two"
                    ),
                    textInput(inputId = "fw_primer_id", 
                        # id for fw primers
                        label = tagList(icon("arrow-right", lib = "glyphicon"),
                        "Identifier for forward primers"), 
                        value = "_fw"
                    ),
                    bsTooltip("fw_primer_id", 
                        paste("The identifier used to declare forward",
                        "primers (complementary to the antisense strand) in the headers of the input FASTA file."),
                        "right", options = list(container = "body")
                    ),
                    textInput(inputId = "rev_primer_id", 
                        # id for reverse primers
                        label = tagList(icon("arrow-left", lib = "glyphicon"),
                        "Identifier for reverse primers"), 
                        value = "_rev"
                    ),
                    bsTooltip("rev_primer_id", 
                        paste("The identifier used to declare reverse primers",
                        "(complementary to the sense strand) in the headers of the input FASTA file."),
                        "right", options = list(container = "body")
                    )#,
                    #####
                    # confirm primers button
                    #######
                    #div(class = "rightAligned",
                        ## header structure confirmation button
                        #actionButton("confirm_primer_header_structure", 
                        #"Confirm primer identifiers", 
                        #icon = icon("check"), 
                        #class="actionStyleSmall", 
                        #style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
                    #)
                ), # close conditional
                conditionalPanel("input.load_eval_primers", 
                    p("There's nothing to configure for loading evaluated primers (CSV input).")
                )
                ) # close panel
                ), # close collapse
                #########
                # primer input:
                #########
                fileInput(inputId = "primer_file", 
                        label = tagList(icon("file"),
                        "Primer FASTA/CSV"), accept="text"
                ),
                bsTooltip("primer_file", 
                    "The set of primers to be evaluated in FASTA format.
                    If primers should be designed from scratch, 
                    no file needs to be uploaded here.",
                    "right", options = list(container = "body")
                )
        ) # close manual primer upload conditional panel
     ), # close conditional panel conditioning on primer evaluation
    conditionalPanel("input.primer_analysis_type == 'design'", 
        ############
        # PRIMER DESIGN PANEL
        ############
        div(p("Please specify the desired properties 
            of the primers that you would like to design.", 
            create.help.button("opti_optimization")),
            class = "two"
        ),
        bsTooltip("help_init_initialization", 
            "View help on initializing the primer set.",
            "right", options = list(container = "body")
        ),
        # tree/naive initialization?
        radioButtons("init_algo", "Template sequence relationship", c("Related" = "tree", "Divergent" = "naive"), inline = TRUE, selected = "tree"),
        bsTooltip("init_algo",
            "If the templates are related a tree-based primer initialization is used, otherwise a naive initialization is employed.",
            "right", options = list(container = "body")
        ),
        # directionality of designed primers
        selectInput("design_direction", 
            "Target strands for design", 
            c("Both strands" = "both", 
              "Antisense strand" = "fw",
              "Sense strand" = "rev")),
        bsTooltip("design_direction",
            "Design forward and reverse primers (\\`Both strands\\`), only forward primers (\\'Antisense\\'), or only reverse primers (\\'Sense\\').",
            "right", options = list(container = "body"),
            trigger = "focus"
        ),

        # degeneracy of primers
        sliderInput("max_degeneracy", "Maximum primer degeneracy", min = 1, max = 64, value = 16),
        bsTooltip("max_degeneracy",
            "Degeneracy is the number of individual oligomers that are represented by the ambiguous sequence of a degenerate primer.",
            "right", options = list(container = "body"),
            trigger = "focus"
        ),

        conditionalPanel("input.init_algo == 'tree'",
            #############
            # TREE INITIALIZATION
            ###########

            # required consevation
            sliderInput("required_conservation", 
                "Percentile of most conserved regions to consider",
                min = 0, max = 1, 
                value = 1
            ),
            bsTooltip("required_conservation", 
                "To improve the runtime, primers can be constructed only in 
                highly-conserved regions. For example, selecting 10% will only
                construct primers in the regions whose conservation is among 
                the top 10% according to Shannon entropy.",
                "right", 
                options = list(container = "body")
            )
        ) 
    ),
    conditionalPanel("input.primer_analysis_type == 'compare'",
        ############
        # COMPARISON PANEL
        ###########
        withTags(
            div(class = "one", 
            p("For immunological applications, we provide a number 
               of evaluated primer sets that can be loaded directly. 
               Otherwise, you can upload raw, evaluated csv files
               containing the primers to be analyzed.", 
            create.help.button("input_primers_comparison"))
        )),
        ##########
        # Choice between supplied/custom comparison primers
        ######
        radioButtons("primer_comparison_upload_choice", 
            tagList(icon("floppy-disk", lib = "glyphicon"), "Primer source"),
            choices = c("Supplied primers" = "available_primers",
                        "Personal primers" = "personal_primers"
                       ), 
            selected = "available_primers", inline=TRUE
        ),
        conditionalPanel("input.primer_comparison_upload_choice == 'personal_primers'",
            ##########
            # PERSONAL PRIMERS
            ############
            # primer csv files for comparison
            fileInput("comparison_file", 
                tagList(icon("tag", lib = "glyphicon"), 
                    "Primer CSV"
                ), multiple = TRUE, accept = "text/csv"
            ),
            bsTooltip("comparison_file", 
                "Upload analyzed, downloaded (raw) primer files 
                 in csv format for comparison.", 
                "right", options = list(container = "body")
            )
            #fileInput("comparison_constraint_files", 
                #tagList(icon("filter", lib = "glyphicon"), 
                    #"Upload constraint settings"
                #), multiple = TRUE, accept = c("text/xml", ".xml")
            #)
        ), # close non-IMGT panel
        conditionalPanel("input.primer_comparison_upload_choice == 'available_primers'",
            #################
            # SUPPLIED COMPARISON PRIMER SETS
            ##################
            # select primer sets for comparison
            selectizeInput("selected_comparison_primers",
            tagList(icon = icon("tag", lib = "glyphicon"),
                    "Available primer sets"),
                 choices = NULL, 
                 selected = NULL, multiple = FALSE, 
                    options = list(
                    placeholder = 'Please select one or multiple of the available primer data sets',
                    onInitialize = I('function() { this.setValue(""); }')
                    )
            ),
            # load comparison primer button
            actionButton("load_all_comparison_sets", 
                "Load all available primer sets",
                class = "actionStyleSmall"
            )
        )
    )
)

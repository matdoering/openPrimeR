###############
# UI panel for inputting templates
##################
tabPanel("Templates", 
    id = "input_data_tabset",
    value = "input_data_panel",
    icon = icon("book", lib = "glyphicon"),
    br(),
    h2(class="inline", "Templates"),
    br(),
    br(),
    # display template introductory text
    div(class = "one", 
        p("Please load a set of template sequences. You can either simply load one of the supplied template data sets or provide your own template data.", 
           create.help.button("input_templates")
        )
    ),
    ### template scenario selector (template source)
    radioButtons("template_scenario", 
                tagList(icon("flask"), "Template source"), 
                choices = c("Available" = "supplied",
                            "Personal" = "personal"), 
                selected = "supplied", inline=TRUE
    ),
    bsTooltip("template_scenario", 
        paste("Load templates provided by openPrimeR or load your own templates?"),
        "right", options = list(container = "body")
    ),
    # selection of type of supplied templates
    conditionalPanel("input.template_scenario == 'supplied'",
         selectizeInput("selected_supplied_templates", 
                        "Template type", 
                        c("Immunological" = "immunological"),
                        selected = "immunological"
         )
    ),
    bsTooltip("selected_supplied_templates", 
        paste("Select one of the available data sets."),
        "right", options = list(container = "body"),
        trigger = "focus"
    ),

    ##############
    # COMPARISON TAB
    ################
    conditionalPanel("input.primer_analysis_type == 'compare'",
        ############
        # template selection for primer comparison
        ###########
        conditionalPanel("input.template_scenario == 'personal'",
            # personal comparison templates
            fileInput("comparison_templates", 
                tagList(icon("file", lib = "glyphicon"), 
                    "Template CSV"
                ), 
                multiple = TRUE, accept = "text/csv"
            ),
            bsTooltip("comparison_templates", 
                      "Upload analyzed, downloaded (raw) template files 
                      in csv format corresponding to the input 
                      primer sets for comparison.", 
                     "right", options = list(container = "body")
            )
        ),
        # IMGT COMPARISON TEMPLATES
        conditionalPanel("input.template_scenario == 'supplied' &&
                         input.selected_supplied_templates == 'immunological'",
            div(class = "one", 
                p("We supply functional template data for IGH,
                   IGK, and IGL of homo sapiens.")
            ), 
            selectizeInput("template_comparison_locus", 
                tagList(icon("book", lib = "glyphicon"), 
                        "Locus"),
                            c("IGH", "IGK", "IGL"),
                            selected = NULL,
                            options = list(
                                placeholder = 'Please select one of the available template data sets',
                                onInitialize = I('function() { this.setValue(""); }')
                            )
            )
        ),
        #######
        # RESET BUTTON FOR COMPARISON TEMPLATES
        #######
        actionButton("reset_rv_comparison.data", "Reset", 
                    icon = icon("refresh", lib = "glyphicon"), 
                    class = "actionStyleSmall"),
        bsTooltip("reset_rv_comparison.data", 
            "Reset uploaded primer and template sets for comparison.", 
            "right", options = list(container = "body")
        )
    ),
    ###########
    # Analysis of primers
    ############
    conditionalPanel("input.primer_analysis_type != 'compare'",
        bsCollapse(id = "template_collapse_analysis",
            open = "template_input_panel",
        bsCollapsePanel(tagList(icon("book", lib = "glyphicon"), 
                        "Template input"),
            value = "template_input_panel", 
            style = "primary",
        ############
        # SUPPLIED TEMPLATES: IMMUNOLOGICAL
        ###############
        conditionalPanel("input.template_scenario == 'supplied' && input.selected_supplied_templates == 'immunological'", 
        # Template upload: other analysis than comparison and IMGT data are used
            div(class = "one", 
                p("Specify the template sequences to be retrieved from ",
                HTML("<a href='http://imgt.org/genedb/' target='_blank'>
                    IMGT Gene-DB</a>"), ". 
                    The allowed primer binding region is automatically set
                    to the leader region."
                )
            ),
            # species to retrieve data for
            selectInput("IMGT_DB_species", 
                "Species", 
                get.IMGT.settings()[["model.gene.id.species.txt"]], 
                selected = "Homo sapiens"),
            # locus to retrieve data for
            selectInput("IMGT_DB_locus", 
                "Locus", 
                get.IMGT.settings()[["model.locusLike.txt"]], 
                selected = "IGH"
            ),
            # IMGT function to retrieve data for
            selectInput("IMGT_DB_function", 
                "Function", 
                get.IMGT.settings()[["model.allele.fcode.txt"]], 
                selected = "functional"
            ),
            # remove partial seqs?
            checkboxInput("remove_partial_seqs", 
                "Remove partial sequences", 
                TRUE
            ), 
            bsTooltip("remove_partial_seqs", 
                "Exclude sequences that are incomplete and augment IMGT data with new sequencing results (if available).",
                "right", options = list(container = "body")
            ),
            # re-load IMGT data stored on disk?
            checkboxInput("update_IMGT_DB_data", 
                "Update existing data", 
                FALSE
            ),
            bsTooltip("update_IMGT_DB_data",
                "Reload data from IMGT?",
                "right", options = list(container = "body")
            ),
            # template retrieve button
            actionButton("IMGT_template_button",
                "Retrieve templates",
                icon = icon("database"), 
                class="actionStyle btn-primary"
            ),
            bsTooltip("IMGT_template_button", 
                      "Retrieve selected templates from IMGT.",
                      "right", options = list(container = "body")
            ),
            # confirm templates  button
            div(class="rightAligned",
                bsButton("IMGT_template_confirm_button",
                    "Confirm templates", 
                    icon = icon("check"), 
                    class="actionStyleSmall", 
                    disabled = TRUE)
            )
            ), # conditional panel for IMGT template ends
            ###############
            # PERSONAL TEMPLATES
            #############
            conditionalPanel("input.template_scenario == 'personal'",
                    ########
                    # TEMPLATE INPUT OPTIONS
                    #########
                    #value = "config_template_header_structure",
                    #style = "primary",
                    div(p("Please customize the following settings according to the data contained in your FASTA file headers and then upload the FASTA file you want to analyze."),
                        class = "two"
                    ),
                    bsCollapse(id = "personal_template_options",
                        open = NULL,
                    #####
                    # start of basic template options
                    ######
                     bsCollapsePanel(
                        tagList(
                            icon("menu-hamburger", lib = "glyphicon"), 
                            "Basic options"),
                        value = "basic_personal_template_options",
                    # header structure choice
                    uiOutput("header_structure", width = "50%"),
                    # header delimiter symbol
                    textInput(inputId = "template_header_delim",
                          label = tagList(icon("list-alt", lib = "glyphicon"), 
                          "Header field delimiter"), 
                          value = "|"
                    ),
                    bsTooltip("header_structure", 
                        "The metadata fields appearing in the headers of the template FASTA file.",
                              "right", 
                        options = list(container = "body"),
                        trigger = "focus"
                    ),
                    bsTooltip("template_header_delim", 
                        "The character separating individual fields in the template headers.",
                              "right", options = list(container = "body")
                    )
                    ),
                    ######
                    # END OF basic options
                    ###########
                  bsCollapsePanel(
                        tagList(
                            icon("menu-hamburger", lib = "glyphicon"), 
                            "Expert options"),
                        value = "template_expert_panel",
                     # column to be used as template identifier
                    selectInput("template_header_ID_column", 
                        tagList(icon("key"), "Header ID field"), 
                        choices  = c("ACCESSION", "GROUP", "SPECIES", "FUNCTION"), 
                        selected = "ACCESSION"
                    ),
                    bsTooltip("template_header_ID_column", "The field in the template header to be used as the identifier of the templates.",
                              "right", options = list(container = "body",
                              trigger = "focus")
                    ),
                    # the character indicating gaps in the templates
                    textInput("gap_char",
                        tagList(icon("scissors"), "Alignment gap character"),
                        value = "-"
                    ),
                    bsTooltip("gap_char", "The character used to indicate gaps in case of aligned input.",
                              "right", options = list(container = "body")
                    ),
                     # remove duplicate seqs?
                    checkboxInput("remove_duplicated_seqs", 
                        "Remove duplicate sequences", 
                        FALSE
                    ), 
                    bsTooltip("remove_duplicated_seqs", 
                        "Exclude sequences that are duplicated.",
                        "right", options = list(container = "body")
                    )
                    ) # end of expert options
                    ), # close options bscollapse/panel
                    br(),
                    bsTooltip("help_input_templates_header", 
                            "View help on defining the template input settings.",
                            "right", options = list(container = "body")
                    ),
                    # PESRSONAL OPTIONS END
                    ##########
                    # START PERSONAL TEMPLATE INPUT
                    ##########
                    # template input panel
                    fileInput(inputId = "sequence_file",
                        label = tagList(icon("file"),
                            "Template FASTA/CSV file"
                        )
                    ),
                    bsTooltip("sequence_file", "The PCR template sequences in FASTA or CSV format.", 
                              "right", options = list(container = "body")
                    ),
                    # confirm templates button
                    div(class = "rightAligned", 
                        bsButton("confirm_uploaded_templates", 
                            "Confirm templates", icon = icon("check"), 
                            class = "actionStyleSmall", disabled = TRUE, 
                            style = "primary"
                        )
                    )
            ) # close conditional panel for personal input
          ), # close panel for template input
############ INSERT START
           bsCollapsePanel(tagList(icon("bookmark"), "Allowed regions"),
                    ########
                    # allowed regions panel
                    #######
                    value = "allowed_template_panel",
                    radioButtons("selected_allowed_region_definition", 
                        label = tagList("Definition of allowed binding region"),
                        choices = c("Template-specific" = "Template-specific" , 
                                    "Uniform" = "Uniform"
                        ), 
                        inline = TRUE 
                    ),
                    bsTooltip("selected_allowed_region_definition", 
                        paste("Define the binding region either individually",
                        "for each template or provide uniform binding intervals"),
                        "right", options = list(container = "body")
                    ),

                    conditionalPanel("input.selected_allowed_region_definition == 'Template-specific'",
                        # Template-specific binding regions
                        div(p("To restrict the allowed sites for primer binding for each template individually, 
                               please input a FASTA file specifying the allowed binding regions for each template.",
                        create.help.button("input_templates_allowed")),
                               class = "two"
                       ),
                       fileInput(inputId = "leader_file",  # input of fw allowed regions file
                            label = tagList(icon("arrow-right", lib = "glyphicon"),
                            "Allowed regions for forward primers (FASTA)")
                        ),
                        # tooltips for fileinputs don't seem to work
                        bsTooltip("leader_file", 
                            paste("A FASTA file with the template binding regions",
                            "(5\\' to 3\\' for forward primers."),
                            "right", options = list(container = "body")
                        ),
                        fileInput(inputId = "leader_file_rev",  # input of rev allowed regions file
                              label = tagList(icon("arrow-left", lib = "glyphicon"),
                              "Allowed regions for reverse primers (FASTA)")
                        ),
                        bsTooltip("leader_file_rev", 
                            paste("A FASTA file with the template binding regions",
                            "(5\\' to 3\\' for reverse primers."),
                                "right", options = list(container = "body")
                        ),
                        ##########
                        ### CONFIRM REGIONS BUTTON:
                        ###########
                        bsTooltip("help_input_templates_allowed", 
                                  "Help on setting the allowed regions for primer binding.",
                                  "right", options = list(container = "body")
                        )
                    ########
                    ),
                    conditionalPanel("input.selected_allowed_region_definition == 'Uniform'",
                        # Uniform binding regions
                        div(p("To restrict the allowed sites for primer binding, please enter the positional range in the templates where the primers should bind. ",
                            create.help.button("input_templates_uniform")),
                            class = "two"
                        ),
                        # allowed regions fw
                        sliderInput("uniform_allowed_regions_fw", 
                            label = tagList(icon("arrow-right", lib = "glyphicon"), 
                            "5' Binding region: forward primers"), 
                            min = 0, 
                            max = 1000, 
                            value = c(1,30)
                        ),
                        bsTooltip("uniform_allowed_regions_fw", 
                            paste("The positional range from the 5\\'",
                                "end of the templates where the forward",
                                "primers should bind."),
                            "right", options = list(container = "body")
                        ),
                        sliderInput("uniform_allowed_regions_rev", 
                            label = tagList(icon("arrow-left", lib = "glyphicon"), 
                            "3' Binding region: reverse primers"), 
                            min = 0,  # if set to 0 -> primer starts with 1
                            max = 1000, 
                            value = c(1,30)
                        ),
                        bsTooltip("uniform_allowed_regions_rev", 
                        paste("The positional range from the 3\\'",
                                "end of the templates where the reverse",
                                "primers should bind."),
                        "right", options = list(container = "body")
                        ),
                        div(class="leftAligned",
                            actionButton("uniform_region_confirm_button", 
                                "Update binding regions", icon = icon("check"), 
                                class="actionStyle btn-primary", 
                            )
                        )
                     ), # uniform panel ends
######################### INSERT HERE##############
 ############
    # CUSTOMIZATION OF REGIONS
    # -> only for template-specific / supplied templates TODO make this part of allowed regions tab
    ########
    conditionalPanel("input.selected_allowed_region_definition == 'Template-specific'",
        # individualization of regions
        #bsCollapse(id = "customize_allowed_regions_collapse", 
            #bsCollapsePanel(tagList(icon("bookmark"), "Customize allowed regions"),
                #value = "customize_allowed_regions_panel",
                div(p("Define the binding range of forward/reverse 
                       primers in relation to the target region."), 
                       class = "two"
                ),
                conditionalPanel("input.individual_allowed_regions_fw[0] != -0.99", # if leaders were loaded correctly
                    # modify fw binding region relative to defined binding region
                    sliderInput("individual_allowed_regions_fw", 
                        label = tagList(icon("arrow-right", lib = "glyphicon"),
                            "Modify the allowed region for forward primers"), 
                        min = -1, 
                        max = 40, 
                        value = c(-0.99, -0.99), step = 1
                    )
                ), 
                bsTooltip("individual_allowed_regions_fw", 
                           paste("The upstream (5\\') binding region for forward primers.",
                           "<br>Negative positions correspond to positions",
                           "upstream of the target amplification site."),
                            "right", options = list(container = "body")
                ),
                conditionalPanel("input.individual_allowed_regions_rev[0] != -0.99",
                    # modify rev binding region relative to defined binding region
                    sliderInput("individual_allowed_regions_rev", 
                        label = tagList(icon("arrow-left", lib = "glyphicon"), "Modify the allowed region for reverse primers"), 
                        min = -1,  # if set to 0 -> primer starts with 1
                        max = 40, value = c(-0.99,-0.99), 
                        step = 1
                    )
                ),
                bsTooltip("individual_allowed_regions_rev", 
                    paste("The downstream (3\\') binding region for reverse primers.",
                           "<br>Negative positions correspond to positions",
                           "downstream of the target amplification site."),
                            "right", options = list(container = "body")
                ),
                # confirm button to delay updates from slides until button is pressed
                div(class="leftAligned",
                    actionButton("individual_region_confirm_button", 
                        "Update binding regions", icon = icon("check"), 
                        class="actionStyle btn-primary"
                    )
                )
            #) # customize region panel ends
            ) # customize conditional ends
            # confirm allowed regions button
            #br()
            #div(class = "rightAligned",
                #bsButton("confirm_uploaded_allowed_regions", 
                        #"Confirm allowed regions", 
                        #icon = icon("check"), 
                        #class = "actionStyleSmall", 
                        #disabled = TRUE, 
                        #style = "primary"
                #)
            #)
          # confirm button end
          ) # allowed regions panel ends
      ) # template collapse ends
    ) # analysis conditional ends
) # template tab ends

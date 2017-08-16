###########
# USER INTERFACE FOR SHINY
###############
# load app dependencies here
library(shiny) # essential
library(shinyBS) # need to attach shinyBS in order to use some functionality..
# load extra functions
source(file.path(system.file("shiny", package = "openPrimeR"), "shiny_server", "extra_IO_shiny.R"))
source(file.path(system.file("shiny", package = "openPrimeR"), "shiny_server", "extra_set_paths.R"))
source(file.path(system.file("shiny", package = "openPrimeR"), "shiny_server", "extra_helper_functions.R"))

ui <- fluidPage(
    shinyjs::useShinyjs(),
    shinyjs::extendShinyjs(text = "shinyjs.reset = function() {history.go(0)}", functions = "reset"),

    includeCSS(file.path(www.folder, "style.css")), # load CSS for help page
    singleton(tags$head(tags$script(src = 'custom_slider.js'))),
    singleton(tags$head(tags$script(src="debounce.js"))),
    # increase shinyBS tooltip delays globally:
    singleton(tags$head(tags$script(src="tooltip_customization.js"))),

    # add layover for progress bar: show when shiny is busy
    singleton(div(
            id = "progressContainer", 
            div(icon("spinner", "fa-spin fa-pulse fa-fw fa-2x"),
                p("Please be patient, computations are running ..."), 
                id = "progressText"
            ),
            div(id = "progressBlock")
        )
    ),
    # reset file input after reload
    singleton(tags$script('
        Shiny.addCustomMessageHandler("resetFileInputHandler", function(x) {   
          var el = $("#" + x);
          el.replaceWith(el = el.clone(true));
          var id = "#" + x + "_progress";     
          $(id).css("visibility", "hidden");
        });
      ')),
        # add JS event to busy modal to re-open it, if we're still busy but the user closes it
    singleton(tags$script('
        $( document ).ready(function() {
        $("#BusyInfo").on("hidden.bs.modal", function (event) {
        if ($("#BusyInfo").hasClass("busy")) {
            $("#BusyInfo").modal("show")
        }
        });
    })
    ')),
    # reduce the size of glyphicons for help
    singleton(tags$head(tags$style(".glyphicon-info-sign {font-size:115%}"))),
    # style for error messages
    singleton(
        tags$head(tags$script(src = "message-handler.js"))
    ),
    singleton(tags$head(
        tags$style(HTML(".shiny-output-error-fatal {color: red;}"))
    )),
    singleton(tags$head(
        tags$style(HTML(".shiny-output-error-critical {color: red;}"))
    )),
    #########
    # STATIC USER MESSAGES
    ##########
    # Tools
    ##########
    bsModal("MissingTools", "Third-party tools", "",
        #htmlOutput("AvailableToolText"),
        DT::dataTableOutput("ToolOverview"),
        size = "large"
    ),
	bsModal("MissingReportDeps", "Report Dependencies Missing", "",
        p("There are missing dependencies (Pandoc/LateX) for creating PDF reports. Please make sure that all required dependencies are installed."),
        size = "small"
    ),
    # Navigation
    ###########
    # design data verification modal:
    singleton(bsModal("DesignVerification", "Design parameter verification", "designButton",
        HTML("<div style='text-align:center;font-size: 14pt'>"), # TODO: space between traffic light? reduce size of checkinput text TODO
        uiOutput("designText"),
        br(),
        traffic_light(),
        br(),
        uiOutput("designTextDiff"),
        HTML("</div>"),
        checkboxInput("evaluate_difficulty_primers", "Estimate number of required primers", value = FALSE),
        actionButton("evaluate_difficulty", "Evaluate problem difficulty", 
            icon = icon("cogs"),
            class="actionStyle btn-primary"),
        HTML("<div style='text-align:center'>"), # TODO: space between traffic light? reduce size of checkinput text TODO
        actionButton("optimizeButton", "Go!", 
            icon = icon("check"),
            class="actionStyleRun btn-primary"),
        HTML("</div>"),
        size = "large")
    ),

    singleton(bsModal("ResetInfo", "Reset session", "",tags$div(icon("refresh", "fa-2x"), 
        p("Do you really want to reset the current session? 
           Click on Reset to reset the app or click on Close to keep the current session.", 
           style = "color:black"
        ), 
        actionButton("reset_button", "Reset", class = "btn-danger"), align = "center"),
        size = "small")
    ),
    singleton(bsModal("ExitInfo", "Exit", "",tags$div(icon("power-off", "fa-2x"), 
                p("Do you really want to exit the app? Click on Exit to exit the app or on Close to remain in the app.", 
                    style = "color:black"
                ), 
                actionButton("exitButton", "Exit", class = "btn-danger"), align = "center"), 
                size = "small"
             )
    ),
    singleton(bsModal("ExitScreen", "Goodbye", "", size = "small", 
                      HTML("<img width = 276, height = 200, src='images/logo.png' alt='logo'/>"), 
                      tags$div(p("Thank you for using openPrimeR. See you again soon!", style = "color:black"), align = "center")
              )
    ),
    #######
    # Data
    #######
    # modal for displaying adapter seqs in primers
    bsModal(id = "AdapterModal", 
              title = "Adapters found", 
              "", # no trigger,
              "The following possible restriction sites were found in the primer set:",
              br(),br(),
               DT::dataTableOutput("primer_restriction_sites"), 
                size = "large"),
    bsModal(id = "NoAdapterModal", 
              title = "No adapters found", 
              "", # no trigger,
              "Your primer set does not seem to contain any restriction sites.",
                size = "large"),
    singleton(bsModal("ProblemEstimationProblem", "Problem difficulty could not be estimated", "",
             tags$p("Could not estimate the problem's difficulty - probably the primer coverage distribution was too narrow.", 
             style = "color:black", align = "center"),
             size = "small")
    ),  # display when user wants to perf
    singleton(bsModal("NotifyNoDataAvailable", "No data available", "",
             tags$p("Could not perform the required action, because no data was available. Please check your input data.", 
             style = "color:black", align = "center"),
             size = "small")
    ),  # display when user wants to perform an action but no data is available
    singleton(bsModal("NotifyIMGT_ConnectionError", "IMGT Data Unavailable", "",
                tags$p("It was not possible to retrieve the selected data from IMGT. 
                        Possible reasons could be having no internet connection or missing dependencies (selenium for python).", 
                        style = "color:black", align = "center"
                ), size = "small")
    ),
    singleton(bsModal("NotifyCouldNotReadFASTA", "Input not readable", "",
              tags$p("Could not read the input file. Please make sure that the file is really in FASTA format.", 
              style = "color:black", align = "center"),
              size = "small")
    ),
    ##########
    # Errors
    ##########
    singleton(bsModal("XML_Parsing_Error", "Could not read settings", "",
              tags$p("Could not read the provided settings XML file. Please check your input!",
                     style = "color:black", align = "center"
              ),
              size = "small")
    ),

    singleton(bsModal("UnexpectedError", "Unexpected Error", "",
              tags$p("An unexpected error occurred. Please help us with fixing the issue by informing us of this 
                     problem by describing when and how the error occured and attaching the error output from the console.", 
                     style = "color:black", align = "center"
              ),
              size = "small")
    ),
    singleton(bsModal("FastaAlphabetError", "Non-supported characters", "",
              tags$p("Some of the input sequences contained non-supported characters. 
                    Please check the console output for more information.", 
                    style = "color:black", align = "center"
              ), 
              size = "small")
    ),
    singleton(bsModal("IDColumnNotFound", "Header ID Column not Found", "",
              tags$p("The specified header ID column could not be found in the header of the templates. Please specify an ID column that is part of the header structure.",
                    style = "color:black", align = "center"
              ), 
              size = "small")
    ),

    ########
    # Templates
    ##########
    singleton(bsModal("TemplateFormatIncorrect", "Incorrect File Format", "",
              tags$p("The structure of the input file did not fulfill the expectations: Please check your input!",
                style = "color:black", align = "center"
              ), 
              size = "small")
    ),
    singleton(bsModal("AllowedRegionTooShort", "Allowed region too short", "",
              tags$p("Could not initialize primers for all template since the allowed binding region was shorter 
                     than the minimal primer length for some templates.", style = "color:black", align = "center"
              ), 
              size = "small")
    ),
    singleton(bsModal("TemplateIDColNotFound", "ID column not found", "",
                      tags$p("The specified ID column was not found in the header of the templates and the first header variable was used as an identifier.", 
                        style = "color:black", align = "center"
                      ),
                      size = "small"
             )
    ),
    singleton(bsModal("TemplateHeaderStructure", "Templates could not be read", "",
                     tags$p("The structure of the header in your template FASTA file did not correspond to the specified structure of the header. 
                            Please check the settings for the header structure specification.", 
                            style = "color:black", align = "center"
                     ), 
                     size = "small"
              )
    ),
    singleton(bsModal("TemplateCoverageUpdateFailed", "Template coverage could not be updated", "",
                     tags$p("The template coverage could not be updated.
                            Please check whether the input primers correspond to the current template sequences.",
                            style = "color:black", align = "center"
                     ), 
                     size = "small"
              )
    ),

    #####
    # Primers
    #####
    singleton(bsModal("NotifyPrimersMissingKeyword", "Primer annotation", "",
                      tags$p("The directionalities of some primers were not explicitly provided. 
                             Some directionalities were inferred and primers could not be paired. Please verify the annotations of the imported primers.", 
                             style = "color:black", align = "center"
                      ), 
                      size = "small"
             )
    ),
    singleton(bsModal("NotifyPrimersNoDirection", "Primer annotation", "",
                     tags$p("Some primers could not be annotated with any directions, hence no primers were imported. 
                            Please check whether you have set the right directionality keywords for your input primers.", 
                            style = "color:black", align = "center"
                     ), 
                     size = "small")
    ),
    singleton(bsModal("NotifyPrimersDuplicateDirections", "Primer annotation", "",
                      tags$p("Some primers had multiple sequences for one direction and therefore no primers were imported. Please check your input files.", 
                        style = "color:black", align = "center"
                      ), 
                      size = "small")
    ),
    #####
    # Binding region
    ####
    singleton(bsModal("AmpliconStartUndefined", "Assignment of allowed regions", "",
                      tags$p("The current assignment of allowed binding regions for primers would not yield any reasonable PCR products for some templates.", 
                      style = "color:black", align = "center"), 
                      size = "small")
    ),
    singleton(bsModal("NotifyAllowedNotFound", "Assignment of allowed regions", "",
                      tags$p("Could not find the provided binding regions in some of the templates. 
                        Please check whether the binding regions agree with the templates.", 
                        style = "color:black", align = "center"
                      ),
                      size = "small")
    ),
    singleton(bsModal("NotifyAllowedNoMatches", "Assignment of allowed regions", "",
                      tags$p("No allowed binding regions could be assigned. 
                             Please check whether the headers of the allowed regions correspond to the template headers.", 
                             style = "color:black", align = "center"
                      ), 
                      size = "small")
     ),
     singleton(bsModal("NotifyAllowedMissing", "Assignment of allowed regions", "",
                        tags$p("Allowed binding regions were not specified for all of the templates
                                and the binding regions of these templates were not adjusted.", 
                                style = "color:black", align = "center"
                       ), 
                       size = "small")
     ),
     singleton(bsModal("NotifyAllowedRedundant", "Assignment of allowed regions", "",
                       tags$p("For some templates, multiple sites for primer binding were specified. 
                            Only the first specified binding region was used for these templates.", 
                            style = "color:black", 
                            align = "center"
                        ),
                        size = "small")
    ),
    singleton(bsModal("NotifyAllowedNotAllLeadersMatched", "Assignment of allowed regions", "",
                      tags$p("At least one of the allowed regions could not be assigned to any template sequence.", 
                        style = "color:black", align = "center"
                      ),
                      size = "small")
    ),
    singleton(bsModal("NotifyNotAllowedBinding", "Binding ratio exceeded", "",
                      tags$p("Warning: The ratio of primers binding to other regions than the allowed binding regions exceeded the allowed ratio.", 
                        style = "color:black", align = "center"
                      ), 
                      size = "small")
    ),
    ##########
    # Optimization
    ####
    singleton(bsModal("RelaxInfoOpti", "Constraint Relaxation", "",
                      tags$p("Constraints were relaxed during the optimization , because target coverage could not be reached with the input constraints.
                          Click on the constraints tab to see the modified constraints.", 
                          style = "color:black", align = "center"
                      ), 
                      size = "small")
    ),

    #######
	# start of UI elements
	#######
    #####
    # header panel: show tool name and image
	#######
    myHeaderPanel(
                tagList(
                        column(4, 
                            div(class="headerStyle", 
                            HTML(paste0("<a href='https://github.molgen.mpg.de/mdoering/primer_design/'", 
                            " target='_blank'><img class='header_img' src='images/logo_text.png' alt='logo'/></a>"))
                        )),
                        column(6, 
                            div(class = "headerStyle", id = "containerHeaderCenter", style = "display:inline-block;",
                                div(
                            selectInput("primer_analysis_type", 
                                tagList(icon = icon("cogs"),
                                    "Analysis mode"
                                ), 
                                c("Evaluation" = "evaluate", 
                                  "Design" = "design", 
                                  "Comparison" = "compare"
                                ),
                                selected = "evaluate",
                                width = "200px"
                            )),
                            div(
                            selectInput("set_meta_selector", 
                                tagList(icon("eye"),
                                        "Set selector"),
                                choices = c("All data" = "all", 
                                            "Filtered data" = "filtered", 
                                            "Design data" = "optimized"),
                                width = "200px"
                            ))
                    )),
                        column(2,
                            div(class="headerStyle", id = "containerHeaderRight",
            
                                div(actionButton("third_party_tools", "", icon = icon("industry"),
                                    class = "actionStyle")),
                                div(actionButton("reset_tool", "", icon = icon("refresh"),
                                    class = "actionStyle")),
                                div(actionButton("quit_tool", "", icon = icon("power-off"),
                                    class = "actionStyle"))
                    )
                )),
                windowTitle = "openPrimeR"
    ),
    br(),
	######	
    # a) show loading message on startup
	######
    div(
        id = "loadingContent",
        div(icon("spinner", "fa-spin fa-pulse fa-fw fa-2x"),
                p("Loading openPrimeR, please be patient ..."), 
                id = "progressTextStartUp"
            )
    ),
	##########
	# b) show true app content
	#########
    shinyjs::hidden(div(id = "app-content",
        ####
        # SIDE PANEL
        ####
        fluidRow(
            column(4,
            wellPanel( # well panel for the sidebar (grey background)
                ##############
                # TEMPLATES TAB
                ################
                tabsetPanel(type = "tabs", id = "settingsPanel",
                    source(file.path(src.ui.folder, "ui_side_templates.R"))$value,
                    ####
                    # INFO SECTION: output important status updates here.
                    ######
                    #br(),
                    #div(textOutput("run_mode", inline=TRUE), style = "color:#214f78;font-weight: bold;"),
                    #div(textOutput("used_nbr_cores", inline=TRUE), style = "color:#214f78;font-weight: bold;")
                    ######
                    # PRIMER TAB
                    #######
                    source(file.path(src.ui.folder, "ui_side_primers.R"))$value,
                    ######
                    # SETTINGS TAB
                    ###### 
                    source(file.path(src.ui.folder, "ui_side_settings.R"))$value,
                    ########
                    # ANALYZE TAB
                    #########
                    source(file.path(src.ui.folder, "ui_side_analyze.R"))$value,
                    #######
                    # DOWNLOAD TAB
                    ########
                    source(file.path(src.ui.folder, "ui_side_download.R"))$value
                    ############
                    # ADD NEW PANELS HERE (if required)
                    #########
            ), # end setings sidebar tabset
            ########
            # NAVIGATION BAR
            ########
            br(),
            br(),
            div(class="rightAligned",
                actionButton("prevBtn", "Previous", icon = icon("arrow-circle-left"), 
                class="actionStyleSmall", style="color: #fff; background-color: #337ab7; border-color: #2e6da4"),
                actionButton("nextBtn", "Next", icon = icon("arrow-circle-right"), 
                class="actionStyleSmall", style="color: #fff; background-color: #337ab7; border-color: #2e6da4")
            ),
            # extend wellPanel with two breaks
            br(),
            br() 
            #####
        ) # end well panel around the sidebar
	), # end column
    ####
    ## START MAIN PANEL
    ####
    column(8,
        tabsetPanel(type = "tabs",
            id = "main",
            ######
            # TEMPLATE DATA
            ######
            source(file.path(src.ui.folder, "ui_main_templates.R"))$value,
            #####
            # PRIMER DATA
            #########
            source(file.path(src.ui.folder, "ui_main_primers.R"))$value,
            #########
            ## COVERAGE DATA
            ########
            source(file.path(src.ui.folder, "ui_main_coverage.R"))$value,
            ###########
            # CONSTRAINT DATA
            ############
            source(file.path(src.ui.folder, "ui_main_constraints.R"))$value,
            ###############
            # COMPARISON DATA
            ###############
            source(file.path(src.ui.folder, "ui_main_comparison.R"))$value,
            ################
            # LOADED SETTINGS
            ################
            source(file.path(src.ui.folder, "ui_main_settings.R"))$value,
            ##########
            # HELP PAGES
            ############
            source(file.path(src.ui.folder, "ui_main_help.R"))$value
        ) # main tabset ends
    ) # main column ends
) # end row
) # end div
) # end hidden
) # end ui

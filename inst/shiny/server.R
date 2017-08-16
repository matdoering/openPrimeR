##########
# Shiny server functionalities
##########

server <- function(input, output, session) {
    source(file.path(system.file("shiny", package = "openPrimeR"), "shiny_server/extra_shiny_backend.R")) # load dependencies for the backend logic of the tool
    shinyjs::hide(selector = "#light") # don't show traffic light for design difficulty when difficulty hasn't been evaluated yet.

    shinyjs::hide(id = "loadingContent", anim = TRUE, animType = "fade") # after dependencies have loaded, hide the loading message
    shinyjs::show("app-content") # show the true app content
    ############
    # convention: reactiveValues (rv) should start with the prefix rv_
    #############
    # rv_values: other general reactive values that do not fit into existing reactive values
    #   relax_info: bsmodal code when filtering relaxation occurred
    #   last_filtering_constraints: last applied filtering constraints
    rv_values <- reactiveValues(
                              "relax_info" = NULL,  
                              "last_filtering_constraints" = NULL 
    )
    ###########################
    # rv_cur input data:
    ###########################
    #   templates_exon: template sequence file
    #   templates_leader: allowed binding regions fw file
    #   templates_leader_rev: allowed binding regions rev file
    #   primers: file with primer sequences
    #   settings: xml file for constraint settings
    rv_cur.input.data <- reactiveValues("templates_exon" = NULL, 
                                        "templates_leader" = NULL,
                                        "templates_leader_rev" = NULL,
                                        "primers" = NULL,
                                        "settings" = NULL) 
    # load all server source files:
    sources <- list.files(system.file("shiny", "shiny_server", 
        package = "openPrimeR"),
        pattern="server_.*.R",
        full.names = TRUE)
    for (s in sources) {
        #message("Loading shiny server source: ", s)
        source(s, local = TRUE)
    }
}
 

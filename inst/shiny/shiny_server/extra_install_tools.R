# Checks and installs all required stand-alone tools
if (!exists("SHINY.PATH")) {
    SHINY.PATH <- system.file("shiny", package = "openPrimeR")
}
source(file.path(SHINY.PATH, "shiny_server", "extra_set_paths.R")) # load paths for tool installation
source(file.path(SHINY.PATH, "shiny_server", "extra_install_helper.R")) # install functions
if (exists("AVAILABLE.TOOLS") && length(AVAILABLE.TOOLS) != 0) {
    AVAILABLE.TOOLS <- install.tools(AVAILABLE.TOOLS)
} else {
    AVAILABLE.TOOLS <- install.tools()
}


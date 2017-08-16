# check dependencies: done here now because docker doesn't use the start script,
# but the ui/server files directly..
###############
server.src.folder <- system.file("shiny", "shiny_server", 
                package = "openPrimeR")
source(file.path(server.src.folder, "extra_set_paths.R"))  # set paths first
#source(file.path(server.src.folder, "extra_install_helper.R")) 
source(file.path(server.src.folder, "extra_IO_shiny.R"))
AVAILABLE.TOOLS <- openPrimeR:::check.tool.function(frontend = TRUE)
source(file.path(server.src.folder, "extra_helper_functions.R"))
# increase font size for ggplot objects for shiny:
OLD.GG.THEME <- ggplot2::theme_get()
ggplot2::theme_set(ggplot2::theme_grey(base_size = 20)) 

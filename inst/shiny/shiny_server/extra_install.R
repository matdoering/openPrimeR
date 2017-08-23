##################
# Checks and installs all R library dependencies of the tool
############################
source(file.path(system.file("shiny", package = "openPrimeR"), "shiny_server", "extra_install_helper.R"))
############ 
cat("###\n# Initializing R libraries ...\n###\n")
check.libPaths()
update.pkgs <- FALSE  # update existing packages? set to FALSE because of shinyjs V8 dependeny on my local machine
if (exists("USED.MIRROR") && length(USED.MIRROR) != 0) {
    # has used.mirror been set in a previous session?
    cat("o Using previously selected repository for package installation ...\n")
    # install.required.packages(USED.MIRROR) # don't need to do this determine
    # USED.MIRROR for later sessions
} else {
    cat("o Fastest repository for package installation will be determined automatically ...\n")
    USED.MIRROR <- install.required.packages(update.pkgs = update.pkgs)
}
# update.packages()

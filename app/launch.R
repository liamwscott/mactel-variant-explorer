# =============================================================================
# launch.R — one-click bootstrap launcher for the MacTel Variant Explorer.
#
# This is the script the double-click launchers call. It:
#   1. Makes the app folder the working directory.
#   2. Checks that every required package is installed; installs any missing
#      ones from CRAN (only happens the first time, or after an R upgrade).
#   3. Starts the Shiny app and opens it in the default web browser.
#
# Collaborators only need R installed (https://cran.r-project.org). They do not
# need RStudio or any command-line knowledge — just double-click the launcher
# for their operating system:
#   macOS    ->  "Run MacTel Explorer.command"
#   Windows  ->  "Run MacTel Explorer.bat"
# =============================================================================

# --- 1. Locate the app folder ------------------------------------------------
# When run via `Rscript launch.R` the launchers already cd into this folder,
# so the working directory is correct. Fall back to the script's own location
# if launched some other way.
get_script_dir <- function() {
  args <- commandArgs(trailingOnly = FALSE)
  file_arg <- grep("^--file=", args, value = TRUE)
  if (length(file_arg) == 1) {
    return(normalizePath(dirname(sub("^--file=", "", file_arg))))
  }
  getwd()
}
app_dir <- get_script_dir()
setwd(app_dir)

# --- 2. Ensure required packages are installed -------------------------------
required <- c("shiny", "bslib", "DT", "ggplot2", "dplyr", "tidyr", "readr",
              "stringr", "forcats", "scales", "plotly", "bsicons", "jsonlite",
              "shinyFiles", "r3dmol", "httr")

installed <- rownames(installed.packages())
missing   <- setdiff(required, installed)

if (length(missing) > 0) {
  message("\n========================================================")
  message("First-time setup: installing ", length(missing), " R package(s):")
  message("  ", paste(missing, collapse = ", "))
  message("This can take a few minutes. It only happens once.")
  message("========================================================\n")
  install.packages(missing, repos = "https://cloud.r-project.org")

  # Re-check; abort with a clear message if anything failed to install.
  still_missing <- setdiff(required, rownames(installed.packages()))
  if (length(still_missing) > 0) {
    stop("Could not install: ", paste(still_missing, collapse = ", "),
         "\nPlease check your internet connection and try again, or install ",
         "these packages manually in R with install.packages().")
  }
}

# --- 3. Launch the app -------------------------------------------------------
message("\nStarting the MacTel Variant Explorer...")
message("It will open in your web browser shortly.")
message("To stop the app: close the browser tab, or close this window, ",
        "or press Esc / Ctrl-C here.\n")

# Quit the R process a couple of seconds after the last browser tab is closed,
# so port 7766 is freed and the next launch starts cleanly (see app.R).
options(mactel.autostop = TRUE)

shiny::runApp(
  appDir        = app_dir,
  port          = 7766,
  host          = "127.0.0.1",
  launch.browser = TRUE
)

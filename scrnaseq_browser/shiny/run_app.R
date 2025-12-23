#!/usr/bin/env Rscript

# Start Shiny app
# Make sure backend is running on http://localhost:8000

message("Starting scRNA-seq Shiny App...")
message("Backend: ", Sys.getenv("API_BASE_URL", "http://localhost:8000"))
message("")

# Check if required packages are installed
required_pkgs <- c(
  "shiny",
  "bslib",
  "httr",
  "jsonlite",
  "base64enc",
  "viridisLite",
  "ggplot2",
  "scales",
  "htmlwidgets",
  "dplyr"
)

missing_pkgs <- required_pkgs[!sapply(required_pkgs, requireNamespace, quietly = TRUE)]
if (length(missing_pkgs) > 0) {
  msg <- paste0(
    "Missing R packages: ", paste(missing_pkgs, collapse = ", "), "\n\n",
    "Install them in your conda env (recommended):\n",
    "  conda install -y -c conda-forge ",
    "r-shiny r-bslib r-httr r-jsonlite r-base64enc r-viridislite r-ggplot2 r-scales r-htmlwidgets r-dplyr\n"
  )
  stop(msg, call. = FALSE)
}

# Respect API_BASE_URL from environment (run_app.sh sets it).
if (Sys.getenv("API_BASE_URL", "") == "") {
  Sys.setenv(API_BASE_URL = "http://localhost:8000")
}

# Respect SHINY_PORT from environment.
shiny_port <- suppressWarnings(as.integer(Sys.getenv("SHINY_PORT", "3838")))
if (is.na(shiny_port) || shiny_port <= 0) shiny_port <- 3838

# Run the app
shiny::runApp(
  appDir = ".",
  host = "0.0.0.0",
  port = shiny_port,
  launch.browser = FALSE
)

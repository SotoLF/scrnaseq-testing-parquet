# Custom scatterplot renderer (no htmlwidgets dependency)

.my_scatterplot_dependency <- function() {
  htmltools::htmlDependency(
    name = "my_scatterplot",
    version = "1.0.0",
    src = c(file = "www"),
    script = "my_scatterplot.js"
  )
}

my_scatterplot <- function(data = NULL, x, y, colorBy = NULL, group_var = NULL,
                           plotId = "scatterplot", size = 3, legend_title = "Group",
                           showAxes = FALSE, width = NULL, height = NULL,
                           elementId = NULL) {

  # Prepare data
  x <- as.numeric(x)
  y <- as.numeric(y)
  points <- data.frame(x = x, y = y)

  # Handle coloring
  color_encoding <- NULL
  if (!is.null(colorBy)) {
    if (is.factor(colorBy) || is.character(colorBy)) {
      # Categorical coloring
      colorBy <- as.factor(colorBy)
      levels_list <- levels(colorBy)
      color_codes <- as.integer(colorBy) - 1
      points$color <- as.integer(color_codes)
      color_encoding <- list(
        type = "categorical",
        levels = levels_list
      )
    } else {
      # Numeric coloring
      points$color <- as.numeric(colorBy)
      color_encoding <- list(
        type = "continuous",
        min = min(colorBy, na.rm = TRUE),
        max = max(colorBy, na.rm = TRUE)
      )
    }
  } else {
    points$color <- NULL
  }

  list(
    plotId = plotId,
    x = points$x,
    y = points$y,
    color = points$color,
    colorEncoding = color_encoding,
    pointSize = size,
    showAxes = showAxes,
    legendTitle = legend_title
  )
}

my_scatterplotOutput <- function(outputId, width = '100%', height = '400px') {
  # NOTE: We intentionally render the widget via renderUI/uiOutput instead of
  # shinyWidgetOutput, to avoid htmlwidgets trying to locate a non-installed
  # package for widget dependencies.
  shiny::div(
    style = sprintf("width:%s; height:%s;", width, height),
    shiny::uiOutput(outputId)
  )
}

renderMy_scatterplot <- function(expr, env = parent.frame(), quoted = FALSE) {
  if (!quoted) { expr <- substitute(expr) }
  force(env)
  shiny::renderUI({
    payload <- eval(expr, envir = env)
    plot_id <- payload$plotId
    if (is.null(plot_id) || !nzchar(plot_id)) plot_id <- "scatterplot"
    el_id <- paste0(plot_id, "_widget")

    dep <- .my_scatterplot_dependency()
    json <- jsonlite::toJSON(payload, auto_unbox = TRUE, digits = 10, null = "null")

    htmltools::tagList(
      dep,
      htmltools::div(id = el_id, style = "width:100%; height:100%;"),
      htmltools::tags$script(htmltools::HTML(sprintf(
        "(function(){\n  var el = document.getElementById(%s);\n  if (!el) return;\n  var payload = %s;\n  if (window.my_scatterplot_render) { window.my_scatterplot_render(el, payload); }\n})();",
        jsonlite::toJSON(el_id, auto_unbox = TRUE),
        json
      )))
    )
  })
}

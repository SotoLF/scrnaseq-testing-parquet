library(shiny)
library(bslib)
library(httr)
library(jsonlite)
library(base64enc)
library(viridisLite)
library(ggplot2)
library(scales)
library(htmlwidgets)
library(reglScatterplot)

# Check for ComplexHeatmap (optional)
has_complex_heatmap <- requireNamespace("ComplexHeatmap", quietly = TRUE)
if (has_complex_heatmap) {
  suppressPackageStartupMessages(library(ComplexHeatmap))
}

options(shiny.fullstacktrace = TRUE)
options(shiny.sanitize.errors = FALSE)

`%||%` <- function(x, y) if (is.null(x)) y else x

API_BASE <- Sys.getenv("API_BASE_URL", "http://localhost:8000")

# ---- encoding helpers (match reglScatterplot payloads) ----
to_base64 <- function(vec) {
  if (is.null(vec)) return(NULL)
  con <- rawConnection(raw(0), "r+")
  writeBin(as.numeric(vec), con, size = 4)
  raw_data <- rawConnectionValue(con)
  close(con)
  paste0("base64:", base64enc::base64encode(raw_data))
}

prepare_gene_payload <- function(vec, name) {
  if (is.null(vec)) return(list(z = NULL, legend = list(var_type = "none")))
  vec <- as.numeric(vec)
  rng <- range(vec, na.rm = TRUE)

  if (!is.finite(rng[1]) || !is.finite(rng[2])) {
    z_norm <- rep(0, length(vec))
    rng <- c(0, 0)
  } else if (rng[2] == rng[1]) {
    z_norm <- rep(0, length(vec))
  } else {
    z_norm <- (vec - rng[1]) / (rng[2] - rng[1])
  }

  cols <- substr(viridisLite::magma(256), 1, 7)
  list(
    z = to_base64(z_norm),
    legend = list(
      minVal = rng[1],
      maxVal = rng[2],
      midVal = mean(vec, na.rm = TRUE),
      var_type = "continuous",
      colors = cols,
      title = name
    )
  )
}

# ---- decoding helpers ----
b64_to_f32 <- function(s) {
  if (is.null(s)) return(NULL)
  if (!startsWith(s, "base64:")) stop("expected base64: prefix")
  raw <- base64enc::base64decode(sub("^base64:", "", s))
  readBin(raw, what="numeric", n=length(raw)/4, size=4, endian="little")
}
b64_to_i32 <- function(s) {
  if (is.null(s)) return(NULL)
  if (!startsWith(s, "base64:")) stop("expected base64: prefix")
  raw <- base64enc::base64decode(sub("^base64:", "", s))
  readBin(raw, what="integer", n=length(raw)/4, size=4, endian="little")
}

# ---- API client ----
api_get <- function(path, query=list(), timeout_sec=30) {
  r <- GET(paste0(API_BASE, path), query=query, timeout(timeout_sec))
  stop_for_status(r)
  fromJSON(content(r, as="text", encoding="UTF-8"), simplifyVector=TRUE)
}
api_post <- function(path, body, timeout_sec=60) {
  r <- POST(paste0(API_BASE, path),
            body=toJSON(body, auto_unbox=TRUE),
            encode="json",
            timeout(timeout_sec))
  stop_for_status(r)
  fromJSON(content(r, as="text", encoding="UTF-8"), simplifyVector=TRUE)
}

FastLoaderAPI <- function(dataset_id) {
  info <- api_get(paste0("/info/", dataset_id))
  list(
    info = function() info,
    datasets = function() api_get("/datasets"),
    scan_obs = function() api_get(paste0("/obs_types/", dataset_id)),
    get_genes = function(measurement="RNA") api_get(paste0("/genes/", dataset_id), query=list(measurement=measurement)),
    get_obsm_xy = function(measurement="RNA", key=NULL) {
      res <- api_get(paste0("/obsm_xy/", dataset_id), query=list(measurement=measurement, key=key))
      data.frame(
        x = b64_to_f32(res$x),
        y = b64_to_f32(res$y)
      )
    },
    get_obs_col = function(col) {
      res <- api_get(paste0("/obs_col/", dataset_id, "/", col))
      if (res$type == "num") {
        return(b64_to_f32(res$data))
      } else {
        codes <- b64_to_i32(res$codes)
        lvls <- res$levels
        indices <- codes + 1
        indices[indices <= 0] <- NA
        return(factor(lvls[indices], levels=lvls))
      }
    },
    get_expression_subset = function(gene, idx0, measurement="RNA", layer="X") {
      res <- api_post("/expr_subset", list(dataset_id=dataset_id, measurement=measurement, layer=layer, gene=gene, idx0=as.integer(idx0)))
      if (is.null(res$data)) return(NULL)
      b64_to_f32(res$data)
    },
    get_expression_batch = function(genes, idx0, measurement="RNA", layer="X") {
      genes_list <- as.list(as.character(genes))
      res <- api_post("/expr_batch", list(dataset_id=dataset_id, measurement=measurement, layer=layer, genes=genes_list, idx0=as.integer(idx0)))
      if (is.null(res$data) || length(res$genes)==0) return(NULL)
      vec <- b64_to_f32(res$data)
      mat <- matrix(vec, nrow=res$n_obs, ncol=res$n_genes, byrow=TRUE)
      colnames(mat) <- res$genes
      mat
    },
    # Optimized: aggregation done in backend
    get_expression_agg = function(genes, group_by, measurement="RNA", layer="X") {
      genes_list <- as.list(as.character(genes))
      res <- api_post("/expr_agg", list(dataset_id=dataset_id, measurement=measurement, layer=layer, genes=genes_list, group_by=group_by), timeout_sec=120)
      if (is.null(res$data) || length(res$genes)==0) return(NULL)
      vec <- b64_to_f32(res$data)
      # Matrix is genes x groups
      mat <- matrix(vec, nrow=res$n_genes, ncol=res$n_groups, byrow=TRUE)
      rownames(mat) <- res$genes
      colnames(mat) <- res$groups
      mat
    }
  )
}

# ---- Scatterplot widget ----
prepare_color_payload <- function(color_vec, title) {
  if (is.null(color_vec) || title == "solid") {
    return(list(colors=NULL))
  }
  if (is.factor(color_vec) || is.character(color_vec)) {
    levels <- levels(as.factor(color_vec))
    if (length(levels) == 0) return(list(colors=NULL))
    cols <- scales::hue_pal()(length(levels))
    names(cols) <- levels
    return(list(colors=cols))
  } else {
    return(list(colors=NULL))
  }
}

# ---- UI ----
ui <- fluidPage(
  tags$head(tags$style(HTML("
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
    body { background: #f8fafc; font-family: 'Inter', sans-serif; color: #1e293b; overflow-x: hidden; }
    .floating-panel { position: fixed; top: 20px; left: 20px; width: 320px; bottom: 20px; background: rgba(255, 255, 255, 0.75);
      backdrop-filter: blur(16px); border: 1px solid rgba(255,255,255,0.8); border-radius: 20px;
      box-shadow: 0 8px 32px rgba(0, 0, 0, 0.05); padding: 24px; overflow-y: auto; z-index: 1000; }
    .main-content { margin-left: 360px; padding: 20px 40px 40px 0; }
    .db-container { max-width: 900px; margin: 80px auto; background: white; padding: 60px; border-radius: 24px; box-shadow: 0 20px 40px rgba(0,0,0,0.03); }
    .db-table { width: 100%; border-collapse: separate; }
    .db-table th { text-align: left; padding: 16px; color: #94a3b8; font-weight: 600; font-size: 12px; border-bottom: 1px solid #e2e8f0; }
    .db-table td { padding: 20px 16px; border-bottom: 1px solid #f8fafc; color: #334155; }
    .btn-load { background: #3b82f6; color: white; border: none; padding: 8px 20px; border-radius: 8px; font-weight: 600; cursor: pointer; }
    .btn-back { margin-bottom: 24px; font-size: 13px; color: #64748b; border: 1px solid #e2e8f0; background: white; padding: 8px 16px; border-radius: 8px; }
    .summary-box { background: linear-gradient(135deg, #6366f1 0%, #8b5cf6 100%); color: white; padding: 24px; border-radius: 16px; text-align: center; margin-bottom: 30px; }
    .big-num { font-size: 36px; font-weight: 800; }
    .plot-wrapper { background: white; border-radius: 12px; box-shadow: 0 4px 6px -1px rgba(0,0,0,0.02); border: 1px solid #f1f5f9; overflow: hidden; height: 100%; position: relative; }
    .nav-underline { display: flex; flex-direction: row; border-bottom: 2px solid #f1f5f9; margin-bottom: 24px; }
    .nav-underline .nav-link { color: #64748b; font-weight: 500; padding: 12px 16px; margin-right: 8px; border-bottom: 2px solid transparent; }
    .nav-underline .nav-link.active { color: #6366f1; border-bottom-color: #6366f1; font-weight: 600; }
    .sidebar-section { margin-bottom: 20px; padding-bottom: 15px; border-bottom: 1px solid #e2e8f0; }
    .sidebar-section:last-child { border-bottom: none; }
    .sidebar-title { font-weight: 600; font-size: 14px; color: #374151; margin-bottom: 12px; }
  "))),
  uiOutput("app_view")
)

server <- function(input, output, session) {

  app_state <- reactiveValues(view="database", dataset_id=NULL)
  z_data <- reactiveValues(obj=NULL, info=NULL, coords=NULL, obs_keys=NULL, loaded_cols=list(), color_map=NULL)

  gene_slots <- reactiveValues(n1=NULL, n2=NULL, n3=NULL)

  datasets_db <- reactiveVal(data.frame(
    id = character(),
    name = character(),
    n_cells = integer(),
    n_genes = integer(),
    stringsAsFactors = FALSE
  ))

  # Load dataset list from API
  observe({
    ds <- tryCatch(
      api_get("/datasets", timeout_sec=10),
      error = function(e) {
        showNotification(paste("Backend error:", conditionMessage(e)), type = "error", duration = 6)
        data.frame(id=character(), name=character(), n_cells=integer(), n_genes=integer(), stringsAsFactors = FALSE)
      }
    )
    datasets_db(ds)
  })

  output$app_view <- renderUI({
    if (app_state$view == "database") {
      div(class="db-container",
          h2("scRNA-seq App", style="font-weight: 800; color: #1e293b;"),
          p("Select a dataset to begin.", style="color:#64748b;"),
          uiOutput("db_list")
      )
    } else {
      req(z_data$info)
      tagList(
        # Floating Sidebar with conditional content based on active tab
        div(class="floating-panel",
            actionButton("go_back", "← Databases", class="btn-back"),
            uiOutput("summary_stats"),
            
            # === VISUALIZATION TAB CONTROLS ===
            conditionalPanel(
              condition = "input.main_tabs == 'Visualization'",
              div(class="sidebar-section",
                  div(class="sidebar-title", "Visualization Settings"),
                  selectInput("p1_color", "Color By", choices=c("solid", z_data$obs_keys), selected="Condition"),
                  selectizeInput("gene_search", "Search Genes (Max 3)", choices=NULL, multiple=TRUE,
                                 options=list(maxItems=3)),
                  sliderInput("pt_size", "Point Size", min=0.5, max=15, value=3, step=0.5),
                  selectInput("layer_select", "Data Layer", choices=z_data$info$layers %||% c("X")),
                  selectInput("emb_select", "Embedding", choices=z_data$info$embeddings %||% c("X_umap"))
              )
            ),
            
            # === VIOLIN TAB CONTROLS ===
            conditionalPanel(
              condition = "input.main_tabs == 'Violin'",
              div(class="sidebar-section",
                  div(class="sidebar-title", "Violin Plot Settings"),
                  selectizeInput("violin_gene", "Gene", choices=NULL),
                  selectInput("violin_group", "Group By", choices=z_data$obs_keys, selected="Condition"),
                  selectInput("violin_layer", "Data Layer", choices=z_data$info$layers %||% c("X")),
                  checkboxInput("violin_points", "Show Points", value=FALSE)
              )
            ),
            
            # === HEATMAP TAB CONTROLS ===
            conditionalPanel(
              condition = "input.main_tabs == 'Heatmap'",
              div(class="sidebar-section",
                  div(class="sidebar-title", "Heatmap Settings"),
                  textAreaInput("heatmap_genes", "Genes (comma separated)", 
                                height = "100px", 
                                placeholder = "CD3D, CD79A, MS4A1, GATA1, HBB..."),
                  selectInput("heatmap_group", "Group By", choices=z_data$obs_keys, selected="Condition"),
                  selectInput("heatmap_layer", "Data Layer", choices=z_data$info$layers %||% c("X")),
                  actionButton("run_heatmap", "Generate Heatmap", class = "btn-load", style="width:100%; margin-top: 10px;")
              )
            )
        ),
        
        # Main Content Area with Tabs
        div(class="main-content",
            navset_underline(
              id="main_tabs",
              
              # TAB 1: VISUALIZATION
              nav_panel("Visualization",
                        fluidRow(
                          column(6, div(class="plot-wrapper", style="height: 40vh; margin-bottom: 20px;",
                                         my_scatterplotOutput("p1", height="100%"))),
                          conditionalPanel(
                            condition = "output.slot1_visible",
                            column(6, div(class="plot-wrapper", style="height: 40vh; margin-bottom: 20px;",
                                           my_scatterplotOutput("p2", height="100%")))
                          ),
                          conditionalPanel(
                            condition = "output.slot2_visible",
                            column(6, div(class="plot-wrapper", style="height: 40vh; margin-bottom: 20px;",
                                           my_scatterplotOutput("p3", height="100%")))
                          ),
                          conditionalPanel(
                            condition = "output.slot3_visible",
                            column(6, div(class="plot-wrapper", style="height: 40vh; margin-bottom: 20px;",
                                           my_scatterplotOutput("p4", height="100%")))
                          )
                        )
              ),
              
              # TAB 2: VIOLIN
              nav_panel("Violin",
                        div(class="plot-wrapper", style="padding: 24px; min-height: 500px;",
                            plotOutput("violin_plot", height="450px")
                        )
              ),
              
              # TAB 3: HEATMAP
              nav_panel("Heatmap",
                        div(class="plot-wrapper", style="padding: 40px; text-align: center; min-height: 600px;",
                            plotOutput("heatmap_plot", height = "550px")
                        )
              )
            )
        )
      )
    }
  })

  output$db_list <- renderUI({
    ds <- datasets_db()
    if (is.null(ds)) {
      return(div(style="color:#64748b; padding: 12px 0;", "Loading datasets…"))
    }
    if (NROW(ds) == 0) {
      return(div(style="color:#64748b; padding: 12px 0;", "No datasets available (backend still starting or unreachable)."))
    }
    rows <- lapply(seq_along(ds$id), function(i) {
      tags$tr(
        tags$td(ds$name[i], style="font-weight:600;"),
        tags$td(ds$id[i]),
        tags$td(tags$button(class="btn-load",
                            onclick=sprintf("Shiny.setInputValue('load_ds', '%s', {priority:'event'})", ds$id[i]),
                            "Analyze"))
      )
    })
    tags$table(class="db-table",
               tags$thead(tags$tr(tags$th("Dataset"), tags$th("ID"), tags$th("Action"))),
               tags$tbody(rows))
  })

  observeEvent(input$load_ds, {
    app_state$dataset_id <- input$load_ds
    showNotification("Connecting to backend API...", type="message", duration=2)

    z <- FastLoaderAPI(app_state$dataset_id)
    info <- z$info()

    default_emb <- if (!is.null(info$embeddings) && length(info$embeddings) > 0) info$embeddings[1] else "X_umap"
    coords <- z$get_obsm_xy(measurement=info$measurement, key=default_emb)
    obs_summary <- z$scan_obs()
    obs_keys <- names(obs_summary)

    z_data$obj <- z
    z_data$info <- info
    z_data$coords <- coords
    z_data$obs_keys <- obs_keys
    z_data$loaded_cols <- list()
    z_data$color_map <- NULL

    gene_slots$n1 <- NULL
    gene_slots$n2 <- NULL
    gene_slots$n3 <- NULL

    app_state$view <- "analysis"
    
    updateSelectInput(session, "layer_select", choices=info$layers, selected="X")
    updateSelectInput(session, "emb_select", choices=info$embeddings, selected=default_emb)
    updateSelectInput(session, "violin_layer", choices=info$layers, selected="X")
    updateSelectInput(session, "heatmap_layer", choices=info$layers, selected="X")
    
    # Find best default for categorical grouping (prefer celltype-like variables)
    celltype_candidates <- c("Annotation", "Cell_Type", "celltype", "cell_type", 
                             "Cluster_Fine", "Cluster_Coarse", "cluster", "leiden", "louvain")
    default_group <- NULL
    for (cand in celltype_candidates) {
      if (cand %in% obs_keys) {
        default_group <- cand
        break
      }
    }
    if (is.null(default_group) && length(obs_keys) > 0) {
      # Fallback: pick first that's not cell_id
      non_id_keys <- obs_keys[!grepl("cell_id|barcode|index", obs_keys, ignore.case = TRUE)]
      default_group <- if (length(non_id_keys) > 0) non_id_keys[1] else obs_keys[1]
    }
    
    updateSelectInput(session, "violin_group", choices=obs_keys, selected=default_group)
    updateSelectInput(session, "heatmap_group", choices=obs_keys, selected=default_group)
    updateSelectInput(session, "p1_color", choices=c("solid", obs_keys), selected=default_group %||% "solid")

    all_genes <- z$get_genes(measurement=info$measurement)
    updateSelectizeInput(session, "gene_search", choices=all_genes, server=TRUE)
    updateSelectizeInput(session, "violin_gene", choices=all_genes, server=TRUE)
    updateSelectizeInput(session, "gene_search", selected=character(0))
  })

  observeEvent(input$emb_select, {
    req(z_data$obj, input$emb_select)
    if (input$emb_select %in% z_data$info$embeddings) {
       coords <- z_data$obj$get_obsm_xy(measurement=z_data$info$measurement, key=input$emb_select)
       z_data$coords <- coords
    }
  })

  observeEvent(input$go_back, {
    app_state$view <- "database"
  })

  ensure_col <- function(col_name) {
    if (is.null(col_name) || col_name == "solid") return(NULL)
    if (is.null(z_data$loaded_cols[[col_name]])) {
      z_data$loaded_cols[[col_name]] <- z_data$obj$get_obs_col(col_name)
    }
    z_data$loaded_cols[[col_name]]
  }

  output$summary_stats <- renderUI({
    req(z_data$info)
    div(class="summary-box",
        div(class="big-num", format(z_data$info$n_cells, big.mark=",")),
        div(class="sub-text", "Total Cells"))
  })

  # ---- gene slots visibility (for conditionalPanel) ----
  .is_valid_gene <- function(g) {
    !is.null(g) && !is.na(g) && nzchar(as.character(g))
  }
  output$slot1_visible <- reactive({ .is_valid_gene(gene_slots$n1) })
  output$slot2_visible <- reactive({ .is_valid_gene(gene_slots$n2) })
  output$slot3_visible <- reactive({ .is_valid_gene(gene_slots$n3) })
  outputOptions(output, "slot1_visible", suspendWhenHidden = FALSE)
  outputOptions(output, "slot2_visible", suspendWhenHidden = FALSE)
  outputOptions(output, "slot3_visible", suspendWhenHidden = FALSE)

  output$p1 <- renderMy_scatterplot({
    req(z_data$coords)
    col_name <- input$p1_color %||% "solid"
    col_vec <- ensure_col(col_name)
    res <- prepare_color_payload(col_vec, col_name)
    if (!is.null(res$colors)) z_data$color_map <- res$colors

    my_scatterplot(
      data=NULL,
      x=z_data$coords$x,
      y=z_data$coords$y,
      colorBy=col_vec,
      group_var=col_vec,
      plotId="p1",
      size=input$pt_size %||% 3,
      legend_title=if (col_name=="solid") "Cells" else col_name,
      showAxes=FALSE
    )
  }, env = environment())

  # Keep point size in sync across plots
  observeEvent(input$pt_size, {
    sz <- input$pt_size %||% 3
    updateMyScatterplotSize(c("p1", "p2", "p3", "p4"), sz)
  }, ignoreInit = TRUE)

  # Fetch up to 3 genes and render as additional linked UMAPs
  observeEvent(list(input$gene_search, input$layer_select), {
    tryCatch({
      req(z_data$obj, z_data$info, z_data$coords)

      genes <- input$gene_search
      if (is.null(genes) || length(genes) == 0) {
        gene_slots$n1 <- NULL; gene_slots$n2 <- NULL; gene_slots$n3 <- NULL
        session$onFlushed(function() {
          session$sendCustomMessage("update_plot_color", list(plotId = "p2", legend = list(var_type = "none")))
          session$sendCustomMessage("update_plot_color", list(plotId = "p3", legend = list(var_type = "none")))
          session$sendCustomMessage("update_plot_color", list(plotId = "p4", legend = list(var_type = "none")))
        }, once = TRUE)
        return()
      }

      genes <- as.character(genes)
      genes <- genes[!is.na(genes) & nzchar(genes)]
      genes <- genes[seq_len(min(3, length(genes)))]
      if (length(genes) == 0) {
        gene_slots$n1 <- NULL; gene_slots$n2 <- NULL; gene_slots$n3 <- NULL
        return()
      }

      layer <- input$layer_select %||% "X"
      idx0 <- seq.int(0, z_data$info$n_cells - 1)

      mat <- z_data$obj$get_expression_batch(genes, idx0, measurement=z_data$info$measurement, layer=layer)
      if (is.null(mat) || ncol(mat) == 0) {
        gene_slots$n1 <- NULL; gene_slots$n2 <- NULL; gene_slots$n3 <- NULL
        return()
      }

      gene_slots$n1 <- genes[1] %||% NULL
      gene_slots$n2 <- if (length(genes) >= 2) genes[2] else NULL
      gene_slots$n3 <- if (length(genes) >= 3) genes[3] else NULL

      g1 <- if (length(genes) >= 1) genes[[1]] else NULL
      g2 <- if (length(genes) >= 2) genes[[2]] else NULL
      g3 <- if (length(genes) >= 3) genes[[3]] else NULL

      session$onFlushed(function() {
        enableMyScatterplotSync(c("p1", "p2", "p3", "p4"), TRUE)

        send_gene <- function(pid, g) {
          if (is.null(g) || is.na(g) || !nzchar(g)) {
            session$sendCustomMessage("update_plot_color", list(plotId = pid, legend = list(var_type = "none")))
            return()
          }
          if (is.null(colnames(mat)) || !(g %in% colnames(mat))) {
            session$sendCustomMessage("update_plot_color", list(plotId = pid, legend = list(var_type = "none")))
            return()
          }
          payload <- prepare_gene_payload(as.numeric(mat[, g]), g)
          session$sendCustomMessage("update_plot_color", list(plotId = pid, z = payload$z, legend = payload$legend))
        }

        send_gene("p2", g1)
        send_gene("p3", g2)
        send_gene("p4", g3)
      }, once = TRUE)
    }, error = function(e) {
      message("[gene_search] ERROR: ", conditionMessage(e))
      gene_slots$n1 <- NULL; gene_slots$n2 <- NULL; gene_slots$n3 <- NULL
      showNotification(paste("Gene plotting error:", conditionMessage(e)), type = "error", duration = 8)
    })
  }, ignoreNULL = FALSE)

  # Render extra plots as solid; gene colors are pushed via `update_plot_color`.
  output$p2 <- renderMy_scatterplot({
    req(z_data$coords)
    my_scatterplot(
      data = NULL,
      x = z_data$coords$x,
      y = z_data$coords$y,
      colorBy = NULL,
      group_var = ensure_col(input$p1_color %||% "solid"),
      plotId = "p2",
      showAxes = FALSE,
      enableDownload = TRUE,
      size = input$pt_size %||% 3
    )
  }, env = environment())

  output$p3 <- renderMy_scatterplot({
    req(z_data$coords)
    my_scatterplot(
      data = NULL,
      x = z_data$coords$x,
      y = z_data$coords$y,
      colorBy = NULL,
      group_var = ensure_col(input$p1_color %||% "solid"),
      plotId = "p3",
      showAxes = FALSE,
      enableDownload = TRUE,
      size = input$pt_size %||% 3
    )
  }, env = environment())

  output$p4 <- renderMy_scatterplot({
    req(z_data$coords)
    my_scatterplot(
      data = NULL,
      x = z_data$coords$x,
      y = z_data$coords$y,
      colorBy = NULL,
      group_var = ensure_col(input$p1_color %||% "solid"),
      plotId = "p4",
      showAxes = FALSE,
      enableDownload = TRUE,
      size = input$pt_size %||% 3
    )
  }, env = environment())

  # ---- VIOLIN PLOT ----
  output$violin_plot <- renderPlot({
    req(z_data$obj, input$violin_gene, input$violin_group)

    gene <- input$violin_gene
    grp_col <- input$violin_group
    layer <- input$violin_layer %||% "X"
    info <- z_data$info

    idx0 <- seq.int(0, info$n_cells - 1)

    expr <- z_data$obj$get_expression_subset(gene, idx0, measurement=info$measurement, layer=layer)
    grp <- ensure_col(grp_col)

    if (is.null(expr) || length(expr) == 0) {
      return(NULL)
    }
    if (is.null(grp) || length(grp) == 0) {
      return(NULL)
    }
    
    if (length(expr) != length(grp)) {
      n <- min(length(expr), length(grp))
      expr <- expr[1:n]
      grp <- grp[1:n]
    }

    df <- data.frame(Expression=expr, Group=as.factor(grp))

    fill_cols <- if (!is.null(z_data$color_map) && grp_col == (input$p1_color %||% "")) z_data$color_map else {
       lvls <- levels(df$Group)
       if (length(lvls) > 0) hue_pal()(length(lvls)) else NULL
    }
    if (!is.null(fill_cols)) names(fill_cols) <- levels(df$Group)

    p <- ggplot(df, aes(x=Group, y=Expression, fill=Group)) +
      geom_violin(scale="width", trim=TRUE, alpha=0.8) +
      theme_minimal() +
      theme(axis.text.x=element_text(angle=45, hjust=1), legend.position="none") +
      labs(title=paste(gene, "Expression by", grp_col), y="Expression", x=NULL) +
      scale_fill_manual(values=fill_cols)

    if (isTRUE(input$violin_points)) {
      p <- p + geom_jitter(size=0.4, alpha=0.3, width=0.2)
    }
    p
  })

  # ---- HEATMAP ----
  # Store heatmap params to detect changes
  heatmap_params <- reactiveValues(
    genes = NULL,
    group = NULL,
    layer = NULL
  )
  
  heatmap_matrix <- eventReactive(input$run_heatmap, {
    req(z_data$obj, input$heatmap_genes, input$heatmap_group)
    
    # Parse genes from comma-separated input
    raw_genes <- strsplit(input$heatmap_genes, ",")[[1]]
    gene_list <- trimws(raw_genes)
    gene_list <- gene_list[gene_list != ""]
    
    if (length(gene_list) == 0) {
      showNotification("Please enter at least one gene", type = "warning")
      return(NULL)
    }
    
    showNotification("Fetching expression data...", type = "message", duration = 2)
    
    # Fetch batch expression data
    info <- z_data$info
    idx0 <- seq.int(0, info$n_cells - 1)
    layer <- input$heatmap_layer %||% "X"
    
    mat <- z_data$obj$get_expression_batch(gene_list, idx0, measurement=info$measurement, layer=layer)
    
    if (is.null(mat) || ncol(mat) == 0) {
      showNotification("No expression data returned for these genes", type = "error")
      return(NULL)
    }
    
    # Fetch grouping variable
    groups <- ensure_col(input$heatmap_group)
    
    if (is.null(groups)) {
      showNotification("Could not load grouping variable", type = "error")
      return(NULL)
    }
    
    # Aggregate by group (mean expression)
    df <- as.data.frame(mat)
    df$Group <- as.character(groups)
    
    agg <- aggregate(. ~ Group, data=df, FUN=mean)
    rownames(agg) <- agg$Group
    agg$Group <- NULL
    
    # Return transposed (Genes x Groups for heatmap)
    result <- t(as.matrix(agg))
    
    showNotification(paste("Heatmap generated:", nrow(result), "genes x", ncol(result), "groups"), type = "message")
    
    return(result)
  })
  
  output$heatmap_plot <- renderPlot({
    mat <- heatmap_matrix()
    req(mat)
    
    t1 <- Sys.time()
    ht <- ComplexHeatmap::Heatmap(mat, 
            name = "Avg Expr",
            col = viridisLite::viridis(100),
            cluster_rows = TRUE,
            cluster_columns = TRUE,
            row_names_gp = grid::gpar(fontsize = 12),
            column_names_gp = grid::gpar(fontsize = 12),
            column_title = paste("Average Expression by", input$heatmap_group),
            rect_gp = grid::gpar(col = "white", lwd = 1)
    )
    ComplexHeatmap::draw(ht)
    t2 <- Sys.time()
    message("Time to render heatmap: ", round(difftime(t2, t1, units="secs"), 2), "s")
  })

}

shinyApp(ui, server)

library(shiny)
library(htmlwidgets)
library(viridisLite)
library(reticulate)
library(bslib)
library(base64enc)
library(reglScatterplot)
library(ggplot2)
library(reshape2)
# Note: ComplexHeatmap is a Bioconductor package. 
# Install: BiocManager::install("ComplexHeatmap")
if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
  # library(ComplexHeatmap)
  suppressPackageStartupMessages(library(ComplexHeatmap))
}

# --- CONFIG & PYTHON ---
# Adjust this path to your specific environment
Sys.setenv(RETICULATE_PYTHON = "/home/lufesu/miniconda3/envs/shiny_app_env/bin/python")
Sys.getenv("RETICULATE_PYTHON")

`%||%` <- function(x, y) if (is.null(x)) y else x

if(file.exists("zarr_loader.py")) {
  source_python("zarr_loader.py") 
} else {
  warning("zarr_loader.py not found.")
}
if(file.exists("my_scatterplot.R")) {
  source("my_scatterplot.R")
}

# --- HELPER FUNCTIONS ---
to_base64 <- function(vec) {
  if (is.null(vec)) return(NULL)
  con <- rawConnection(raw(0), "r+")
  writeBin(as.numeric(vec), con, size = 4)
  raw_data <- rawConnectionValue(con)
  close(con)
  paste0("base64:", base64enc::base64encode(raw_data))
}

# Now returns BOTH the payload and the named color vector for R-side use
prepare_color_payload <- function(color_vec, title) {
  if (is.null(color_vec) || title == "solid") {
     return(list(
       payload = list(z = NULL, legend = list(var_type="none")),
       colors = NULL
     ))
  }
  
  if (is.character(color_vec) || is.factor(color_vec)) {
      levels <- levels(as.factor(color_vec))
      cols <- if (requireNamespace("scales", quietly = TRUE)) scales::hue_pal()(length(levels)) else rainbow(length(levels))
      names(cols) <- levels
      z_norm <- as.integer(as.factor(color_vec)) - 1L
      
      return(list(
         payload = list(z = to_base64(z_norm), legend = list(names = levels, colors = as.vector(cols), var_type = "categorical", title = title)),
         colors = cols # Named vector
      ))
  } else {
      rng <- range(color_vec, na.rm = TRUE)
      z_norm <- (color_vec - rng[1]) / (rng[2] - rng[1])
      cols <- substr(viridisLite::viridis(256), 1, 7)
      
      return(list(
         payload = list(z = to_base64(z_norm), legend = list(minVal = rng[1], maxVal = rng[2], midVal = mean(rng), var_type = "continuous", colors = cols, title = title)),
         colors = NULL
      ))
  }
}

prepare_gene_payload <- function(vec, name) {
  if (is.null(vec)) return(list(z = NULL, legend = NULL))
  rng <- range(vec, na.rm = TRUE)
  if (rng[2] == rng[1]) {
    z_norm <- rep(0, length(vec))
  } else {
    z_norm <- (vec - rng[1]) / (rng[2] - rng[1])
  }
  cols <- substr(viridisLite::magma(256), 1, 7)
  return(list(
    z = to_base64(z_norm), 
    legend = list(minVal = rng[1], maxVal = rng[2], midVal = mean(rng), var_type = "continuous", colors = cols, title = name)
  ))
}

FastLoader <- function(url) {
  py_obj <- py$ZarrBackend(url)
  return(list(
    info = function() py_obj$list_info(),
    scan_obs = function() py_obj$scan_obs(),
    get_obsm = function(key) {
      mat <- py_obj$get_obsm(key)
      df <- as.data.frame(mat); colnames(df) <- c("x", "y")
      return(df)
    },
    get_obs_col = function(name) {
      res <- py_obj$get_obs(name)
      if(is.null(res)) return(NULL)
      unlist(res)
    },
    get_expression = function(gene, layer="X") {
      res <- py_obj$get_expression(gene, layer)
      if(is.null(res)) return(NULL)
      as.numeric(res)
    },
    # [NEW] Batch Fetcher for Heatmap
    get_expression_batch = function(genes, layer="X") {
      res <- py_obj$get_expression_batch(genes, layer)
      if(is.null(res)) return(NULL)
      # Returns matrix (cells x genes)
      return(res)
    },
    get_genes = function() as.vector(py_obj$genes)
  ))
}

# --- DATASETS DB ---
datasets_db <- data.frame(
  Name = c("Bone Marrow (May 2024)", "Corrected Doublet DecontX", "Strict EPDSC Annotated"),
  URL = c("https://scrnaseq-browser.s3.us-east-2.amazonaws.com/BoneMarrow_May2024_outliers_removed_Palantir.zarr",
          "https://scrnaseq-browser.s3.us-east-2.amazonaws.com/corrected_doublet_decontX_mt_Duo_allergy.zarr",
          "https://scrnaseq-browser.s3.us-east-2.amazonaws.com/strict_epdsc_annotated_data_csc_full.zarr"),
  Cells = c("75k", "~2k", "~7.5k"),
  Genome = c("mm10", "mm10", "mm10"),
  stringsAsFactors = FALSE
)

# --- UI ---
ui <- fluidPage(
  tags$head(tags$style(HTML("
    @import url('https://fonts.googleapis.com/css2?family=Inter:wght@400;500;600;700&display=swap');
    body { background: #f8fafc; font-family: 'Inter', sans-serif; color: #1e293b; overflow-x: hidden; }
    .floating-panel {
      position: fixed; top: 20px; left: 20px; width: 320px; bottom: 20px;
      background: rgba(255, 255, 255, 0.75);
      backdrop-filter: blur(16px);
      border: 1px solid rgba(255,255,255,0.8);
      border-radius: 20px;
      box-shadow: 0 8px 32px rgba(0, 0, 0, 0.05);
      padding: 24px; overflow-y: auto; z-index: 1000;
    }
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
  "))),

  uiOutput("app_view")
)

# --- SERVER ---
server <- function(input, output, session) {
  
  # --- STATE ---
  app_state <- reactiveValues(view = "database", db_index = NULL)
  
  z_data <- reactiveValues(
      obj=NULL, info=NULL, coords=NULL, 
      cat_keys=NULL, num_keys=NULL, loaded_cols=list(),
      color_map = NULL # Global Color Sync
  )
  
  gene_slots <- reactiveValues(n1=NULL, n2=NULL, n3=NULL)
  gene_data <- reactiveValues(d1=NULL, d2=NULL, d3=NULL)
  
  # --- VIEW CONTROLLER ---
  output$app_view <- renderUI({
    if (app_state$view == "database") {
      # === DATABASE VIEW ===
      div(class="db-container",
        h2("scRNA-seq App", style="font-weight: 800; color: #1e293b;"),
        p("Select a dataset to begin.", style="color:#64748b;"),
        uiOutput("db_list")
      )
    } else {
      # === ANALYSIS VIEW ===
      req(z_data$info)
      
      tagList(
        # 1. Floating Sidebar
        div(class = "floating-panel",
            actionButton("go_back", "← Databases", class="btn-back"),
            uiOutput("summary_stats"),
            
            # --- Sidebar Content Changes Based on Tab ---
            conditionalPanel(
              condition = "input.main_tabs == 'Visualization'",
              selectInput("p1_color", "Annotation", choices = c("Solid Color" = "solid", z_data$cat_keys)),
              selectizeInput("gene_search", "Search Genes (Max 3)", choices = NULL, multiple = TRUE, options = list(maxItems = 3)),
              sliderInput("pt_size", "Point Size", min=0.5, max=15, value=3, step=0.5),
              selectInput("layer_select", "Data Layer", choices = if(!is.null(z_data$info$layers)) z_data$info$layers else c("X")),
              hr(),
              selectizeInput("sel_meta", "Add Filter", choices = z_data$num_keys, multiple = TRUE),
              uiOutput("dynamic_filters")
            ),
            
            conditionalPanel(
              condition = "input.main_tabs == 'Violin'",
              h5("Violin Settings", style="font-weight:600;"),
              selectInput("violin_gene", "Gene", choices = NULL), # Filled server-side
              selectInput("violin_group", "Group By", choices = z_data$cat_keys),
              selectInput("violin_layer", "Data Layer", choices = if(!is.null(z_data$info$layers)) z_data$info$layers else c("X")),
              checkboxInput("violin_points", "Show Points", value = FALSE)
            ),
            
            conditionalPanel(
              condition = "input.main_tabs == 'Heatmap'",
              h5("Heatmap Settings", style="font-weight:600;"),
              textAreaInput("heatmap_genes", "Genes (comma separated)", height = "100px", placeholder = "CD3D, CD79A, MS4A1..."),
              selectInput("heatmap_group", "Group By", choices = z_data$cat_keys),
              selectInput("heatmap_layer", "Data Layer", choices = if(!is.null(z_data$info$layers)) z_data$info$layers else c("X")),
              actionButton("run_heatmap", "Generate Heatmap", class = "btn-load", style="width:100%")
            )
        ),
        
        # 2. Main Workspace
        div(class = "main-content",
            navset_underline(
              id = "main_tabs",
              
              # TAB 1: VISUALIZATION
              nav_panel("Visualization",
                fluidRow(
                  column(6, div(class="plot-wrapper", style="height: 40vh; margin-bottom: 20px;", my_scatterplotOutput("p1", height="100%"))),
                  conditionalPanel(condition = "output.slot1_visible", column(6, div(class="plot-wrapper", style="height: 40vh; margin-bottom: 20px;", my_scatterplotOutput("p2", height="100%")))),
                  conditionalPanel(condition = "output.slot2_visible", column(6, div(class="plot-wrapper", style="height: 40vh; margin-bottom: 20px;", my_scatterplotOutput("p3", height="100%")))),
                  conditionalPanel(condition = "output.slot3_visible", column(6, div(class="plot-wrapper", style="height: 40vh; margin-bottom: 20px;", my_scatterplotOutput("p4", height="100%"))))
                )
              ),
              
              # TAB 2: VIOLIN
              nav_panel("Violin",
                 div(class="plot-wrapper", style="padding: 24px; height: 600px;",
                     plotOutput("violin_plot", height = "100%")
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

  # --- OUTPUT VISIBILITY ---
  output$slot1_visible <- reactive({ !is.null(gene_slots$n1) })
  output$slot2_visible <- reactive({ !is.null(gene_slots$n2) })
  output$slot3_visible <- reactive({ !is.null(gene_slots$n3) })
  outputOptions(output, "slot1_visible", suspendWhenHidden = FALSE)
  outputOptions(output, "slot2_visible", suspendWhenHidden = FALSE)
  outputOptions(output, "slot3_visible", suspendWhenHidden = FALSE)

  # --- DB LIST ---
  output$db_list <- renderUI({
    rows <- lapply(1:nrow(datasets_db), function(i) {
      tags$tr(
        tags$td(datasets_db$Name[i], style="font-weight:600;"),
        tags$td(datasets_db$Cells[i]),
        tags$td(datasets_db$Genome[i]),
        tags$td(tags$button(class = "btn-load", onclick = sprintf("Shiny.setInputValue('load_db_trigger', %d, {priority: 'event'})", i), "Analyze"))
      )
    })
    tags$table(class = "db-table", tags$thead(tags$tr(tags$th("Dataset"), tags$th("Cells"), tags$th("Genome"), tags$th("Action"))), tags$tbody(rows))
  })

  # --- LOAD LOGIC ---
  observeEvent(input$load_db_trigger, {
    idx <- as.numeric(input$load_db_trigger)
    app_state$db_index <- idx
    
    showNotification("Connecting to Zarr...", type = "message", duration = 3)
    
    url <- datasets_db$URL[idx]
    z <- FastLoader(url)
    info <- z$info()
    
    umap_key <- grep("umap", info$obsm, value = TRUE, ignore.case = TRUE)[1]
    if(is.na(umap_key)) umap_key <- info$obsm[1]
    coords <- z$get_obsm(umap_key)
    
    obs_summary <- z$scan_obs()
    cat_keys <- names(obs_summary)[obs_summary == "cat"]
    num_keys <- names(obs_summary)[obs_summary == "num"]
    
    n_pts <- if(!is.null(info$n_cells)) info$n_cells else nrow(coords)
    opt_size <- max(1.5, min(6.0, round((600 / sqrt(n_pts)) * 0.8, 1)))
    updateSliderInput(session, "pt_size", value = opt_size)
    
    z_data$obj <- z; z_data$info <- info; z_data$coords <- coords
    z_data$cat_keys <- cat_keys; z_data$num_keys <- num_keys
    z_data$loaded_cols <- list(); z_data$color_map <- NULL
    
    app_state$view <- "analysis"
  })

  observeEvent(input$go_back, { app_state$view <- "database" })
  
  observeEvent(app_state$view, {
    req(app_state$view == "analysis", z_data$obj)
    default_cat <- "solid"
    for(c in c("Cell_Type", "Annotation", "leiden", "louvain")) { if(c %in% z_data$cat_keys) { default_cat <- c; break } }
    
    all_genes <- z_data$obj$get_genes()
    updateSelectInput(session, "p1_color", selected = default_cat)
    updateSelectInput(session, "violin_group", selected = default_cat)
    updateSelectInput(session, "heatmap_group", selected = default_cat)
    updateSelectizeInput(session, "gene_search", choices = all_genes, server = TRUE)
    updateSelectizeInput(session, "violin_gene", choices = all_genes, server = TRUE)
  }, once = FALSE)

  ensure_col <- function(col_name) {
    if (is.null(col_name) || col_name == "solid") return(NULL)
    if (is.null(z_data$loaded_cols[[col_name]])) {
       z_data$loaded_cols[[col_name]] <- z_data$obj$get_obs_col(col_name)
    }
    return(z_data$loaded_cols[[col_name]])
  }

  # --- GLOBAL SUBSET LOGIC ---
  subset_indices <- reactive({
    req(z_data$obj)
    total_cells <- z_data$info$n_cells
    valid <- rep(TRUE, total_cells)
    
    if (!is.null(input$sel_meta)) {
      for (col in input$sel_meta) {
        rng <- input[[paste0("rng_", col)]]
        if (!is.null(rng)) {
          vec <- ensure_col(col)
          if (!is.null(vec)) valid <- valid & (vec >= rng[1] & vec <= rng[2])
        }
      }
    }
    
    # Lasso Selection
    sel_data <- input$p1_selected
    if (!is.null(sel_data) && !is.null(sel_data$indices) && length(sel_data$indices) > 0) {
      sel_indices <- unlist(sel_data$indices) + 1 
      lasso_mask <- rep(FALSE, total_cells)
      lasso_mask[sel_indices] <- TRUE
      valid <- valid & lasso_mask
    }
    
    which(valid)
  })

  output$summary_stats <- renderUI({
    req(z_data$info)
    current_indices <- subset_indices()
    n_subset <- length(current_indices)
    pct <- round((n_subset / z_data$info$n_cells) * 100, 1)
    div(class = "summary-box", div(class="big-num", format(n_subset, big.mark=",")), div(class="sub-text", paste0(pct, "% of Total Cells")))
  })

  # --- MAIN PLOT (Updated for Filtering) ---
  output$p1 <- renderMy_scatterplot({
    req(z_data$coords)
    
    # 1. Listen to the subset changes so the plot updates when filters change
    idx <- subset_indices() 
    
    isolate({
        col_name <- input$p1_color %||% "solid"
        col_vec <- ensure_col(col_name)
        
        # 2. Subset the data before sending to the plot
        # We must subset x, y, and the color vector
        x_sub <- z_data$coords$x[idx]
        y_sub <- z_data$coords$y[idx]
        
        # Handle Color Subsetting
        if (!is.null(col_vec) && length(col_vec) == nrow(z_data$coords)) {
          col_vec_sub <- col_vec[idx]
        } else {
          col_vec_sub <- col_vec # Handle cases where it might be NULL or simple string
        }
        
        # Recalculate color payload for the subset
        res <- prepare_color_payload(col_vec_sub, col_name)
        
        # 3. Update global map if needed (Optional: careful with subsets affecting global legends)
        # if(!is.null(res$colors)) z_data$color_map <- res$colors 

        sz <- if(!is.null(input$pt_size)) input$pt_size else 3
        
        my_scatterplot(
          data = NULL, 
          x = x_sub, 
          y = y_sub, 
          colorBy = col_vec_sub, 
          group_var = col_vec_sub, 
          plotId = "p1", 
          size = sz, 
          legend_title = if(col_name=="solid") "Cells" else col_name, 
          showAxes = FALSE
        )
    })
  })

  observeEvent(input$p1_color, {
    req(z_data$obj)
    col_name <- input$p1_color
    col_vec <- ensure_col(col_name)
    res <- prepare_color_payload(col_vec, col_name) 
    z_data$color_map <- res$colors 
    session$sendCustomMessage("update_plot_color", list(plotId = "p1", z = res$payload$z, group_data = res$payload$z, legend = res$payload$legend))
  }, ignoreInit = TRUE)
  
  observeEvent(input$pt_size, {
     updateMyScatterplotSize(c("p1", "p2", "p3", "p4"), input$pt_size)
  }, ignoreInit = TRUE)

  observeEvent(list(input$gene_search, input$layer_select), {
    req(z_data$obj)
    vec <- input$gene_search; layer <- input$layer_select %||% "X"
    update_slot <- function(i, slot_name, data_name, pid) {
      new_val <- if(length(vec) >= i) vec[i] else NULL
      gene_slots[[slot_name]] <- new_val
      if(is.null(new_val)) {
         gene_data[[data_name]] <- NULL
      } else {
         val <- z_data$obj$get_expression(new_val, layer)
         gene_data[[data_name]] <- val
         payload <- prepare_gene_payload(val, new_val)
         session$sendCustomMessage("update_plot_color", list(plotId = pid, z = payload$z, legend = payload$legend))
      }
    }
    update_slot(1, "n1", "d1", "p2")
    update_slot(2, "n2", "d2", "p3")
    update_slot(3, "n3", "d3", "p4")
    enableMyScatterplotSync(c("p1", "p2", "p3", "p4"), TRUE)
  }, ignoreNULL = FALSE)

  render_gene <- function(name_key, data_key, pid) {
    renderMy_scatterplot({
      g_name <- gene_slots[[name_key]]
      g_data <- gene_data[[data_key]]
      
      if(is.null(g_name) || is.null(g_data)) return(NULL)
      
      # 1. Get Indices
      idx <- subset_indices()
      
      isolate({
          col_vec <- ensure_col(input$p1_color)
          sz <- if(!is.null(input$pt_size)) input$pt_size else 3
          
          # 2. Subset Data
          x_sub <- z_data$coords$x[idx]
          y_sub <- z_data$coords$y[idx]
          g_data_sub <- g_data[idx] # Subset expression data
          
          # Subset grouping var if used for tooltips/etc
          col_vec_sub <- if(!is.null(col_vec)) col_vec[idx] else NULL

          my_scatterplot(
            data = NULL, 
            x = x_sub, 
            y = y_sub, 
            colorBy = g_data_sub, 
            group_var = col_vec_sub, 
            continuous_palette = "magma", 
            plotId = pid, 
            legend_title = g_name, 
            showAxes = FALSE, 
            enableDownload = TRUE, 
            size = sz
          )
      })
    })
  }
  output$p2 <- render_gene("n1", "d1", "p2")
  output$p3 <- render_gene("n2", "d2", "p3")
  output$p4 <- render_gene("n3", "d3", "p4")

  output$dynamic_filters <- renderUI({
    req(input$sel_meta)
    lapply(input$sel_meta, function(v) {
      vec <- ensure_col(v); rng <- round(range(vec, na.rm=TRUE), 2)
      div(style="margin-bottom:15px;",
          div(style="font-weight:600; font-size:11px; text-transform:uppercase; color:#94a3b8;", v),
          plotOutput(paste0("hist_", v), height="40px"),
          sliderInput(paste0("rng_", v), NULL, min=rng[1], max=rng[2], value=rng, ticks=FALSE, width="100%")
      )
    })
  })
  
  observe({
    req(input$sel_meta)
    for(v in input$sel_meta) {
      local({
        my_v <- v; vec <- ensure_col(my_v)
        output[[paste0("hist_", my_v)]] <- renderPlot({
          par(mar=c(0,0,0,0), bg=NA); hist(vec, main=NA, axes=FALSE, col="#cbd5e1", border="white", breaks=20)
        }, bg="transparent")
      })
    }
  })

  # --- VIOLIN PLOT IMPLEMENTATION ---
  output$violin_plot <- renderPlot({
    req(z_data$obj, input$violin_gene, input$violin_group)
    
    idx <- subset_indices()
    if (length(idx) < 5) return(NULL)
    
    gene <- input$violin_gene
    group_col <- input$violin_group
    
    # Fetch Data
    expr <- z_data$obj$get_expression(gene, input$violin_layer)
    groups <- ensure_col(group_col)
    
    if (is.null(expr) || is.null(groups)) return(NULL)
    
    # Subset
    plot_df <- data.frame(
      Expression = expr[idx],
      Group = as.factor(groups[idx])
    )
    
    # Determine colors from Global Map
    fill_cols <- if (!is.null(z_data$color_map) && group_col == input$p1_color) {
      z_data$color_map
    } else {
      # Fallback palette
      scales::hue_pal()(length(levels(plot_df$Group)))
    }
    
    p <- ggplot(plot_df, aes(x=Group, y=Expression, fill=Group)) +
      geom_violin(scale = "width", trim = TRUE, alpha=0.8) +
      theme_minimal() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none") +
      labs(title = paste(gene, "Expression by", group_col), y = "Expression", x = NULL) +
      scale_fill_manual(values = fill_cols)
      
    if(input$violin_points) {
      p <- p + geom_jitter(size=0.5, alpha=0.4, width=0.2)
    }
    
    p
  })

  # --- HEATMAP IMPLEMENTATION ---
  
  # Reactive to fetch heatmap data
  heatmap_matrix <- eventReactive(input$run_heatmap, {
    req(z_data$obj, input$heatmap_genes, input$heatmap_group)
    
    # 1. Parse Genes
    raw_genes <- strsplit(input$heatmap_genes, ",")[[1]]
    gene_list <- trimws(raw_genes)
    gene_list <- gene_list[gene_list != ""]
    
    if(length(gene_list) == 0) return(NULL)
    
    # 2. Fetch Batch Data (Cells x Genes)
    mat <- z_data$obj$get_expression_batch(gene_list, input$heatmap_layer)
    if(is.null(mat)) return(NULL)
    
    # Add column names
    # Note: get_expression_batch filters invalid genes, so we need to know which ones returned.
    # The python function returns valid genes in order of indices.
    # For simplicity here, we assume user typed correctly or we map based on returned shape.
    # A robust app would have Python return the valid gene names too. 
    # Here we assume standard behavior:
    colnames(mat) <- gene_list[1:ncol(mat)] 
    
    # 3. Fetch Grouping
    groups <- ensure_col(input$heatmap_group)
    idx <- subset_indices()
    
    # Subset
    mat_sub <- mat[idx, , drop=FALSE]
    groups_sub <- groups[idx]
    
    # 4. Aggregate (Mean by Group)
    # Aggregating is safer for web apps than showing 50k cells
    df <- as.data.frame(mat_sub)
    df$Group <- groups_sub
    
    agg <- aggregate(. ~ Group, data=df, FUN=mean)
    rownames(agg) <- agg$Group
    agg$Group <- NULL
    
    # Return as Matrix (Genes x Groups for ComplexHeatmap)
    return(t(as.matrix(agg)))
  })

  output$heatmap_plot <- renderPlot({
    mat <- heatmap_matrix()
    req(mat)
    
    # Render with ComplexHeatmap if available, else standard heatmap
    if (requireNamespace("ComplexHeatmap", quietly = TRUE)) {
      ComplexHeatmap::Heatmap(mat, 
              name = "Avg Expr",
              col = viridisLite::viridis(100),
              cluster_rows = TRUE,
              cluster_columns = TRUE,
              row_names_gp = grid::gpar(fontsize = 12),
              column_names_gp = grid::gpar(fontsize = 12),
              column_title = paste("Average Expression by", input$heatmap_group),
              rect_gp = grid::gpar(col = "white", lwd = 1)
      )
    } else {
      # Fallback to base R
      heatmap(mat, scale="row", col=cm.colors(256))
    }
  })
}

shinyApp(ui, server)

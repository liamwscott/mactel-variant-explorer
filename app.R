# =============================================================================
# MacTel Variant Explorer — Shiny app
#
# Interactive exploration of Cavalier candidate-variant output.
# Expects the 71-column Cavalier CSV structure (see R/load_data.R).
#
# Run locally:
#   shiny::runApp("/Users/scott.l/Documents/Claude/mactel_variant_explorer")
# or open app.R in RStudio and click "Run App".
#
# Data source priority:
#   1. data/candidate_variants.csv  (real data, gitignored — used if present)
#   2. data/example_variants.csv    (de-identified example shipped with repo)
#   3. a file uploaded via the sidebar
# =============================================================================

library(shiny)
library(bslib)
library(DT)
library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)
library(stringr)
library(forcats)
library(scales)
library(plotly)

# Files in R/ are auto-sourced by Shiny, but source explicitly so the app also
# works when launched via shiny::runApp() from any working directory.
app_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
for (f in list.files(file.path(app_dir, "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

DEFAULT_REAL    <- file.path(app_dir, "data", "candidate_variants.csv")
DEFAULT_EXAMPLE <- file.path(app_dir, "data", "example_variants.csv")
startup_path    <- if (file.exists(DEFAULT_REAL)) DEFAULT_REAL else DEFAULT_EXAMPLE

# Gene info (Tier + descriptions; gene symbols only, no participant data) ------
GENE_INFO_PATH <- file.path(app_dir, "data", "gene_info.tsv")
GENE_INFO      <- load_gene_info(GENE_INFO_PATH)

# Pfam protein domains, pre-fetched offline (gene symbols only, no patient data).
PROTEIN_DOMAINS <- load_protein_domains(file.path(app_dir, "data",
                                                  "protein_domains.tsv"))

# Tier lookup: prefer the richer gene_info table; fall back to gene_tiers.tsv.
TIER_PATH <- file.path(app_dir, "data", "gene_tiers.tsv")
TIER_DF   <- if (!is.null(GENE_INFO)) {
  dplyr::distinct(dplyr::select(GENE_INFO, SYMBOL, Tier))
} else {
  load_gene_tiers(TIER_PATH)
}

# Build a modalDialog describing a single gene from GENE_INFO.
show_gene_modal <- function(symbol) {
  if (is.null(symbol) || is.na(symbol) || symbol == "") return(invisible())
  info <- if (!is.null(GENE_INFO)) GENE_INFO[GENE_INFO$SYMBOL == symbol, ] else NULL

  if (is.null(info) || nrow(info) == 0) {
    body <- tags$p(tags$em("No annotation available for this gene."))
    tier <- NULL
  } else {
    info <- info[1, ]
    fld <- function(label, value) {
      if (is.null(value) || is.na(value) || value == "") return(NULL)
      tags$p(tags$strong(paste0(label, ": ")), value)
    }
    body <- tagList(
      fld("Ensembl ID",       info$Ensembl_ID),
      fld("Chromosome",       info$Chromosome),
      fld("Evidence category", info$Evidence_Category),
      fld("Evidence detail",  info$Evidence_Detail),
      if (!is.null(info$Gene_Description) && !is.na(info$Gene_Description) &&
          info$Gene_Description != "") {
        tagList(tags$hr(),
                tags$p(tags$strong("Description")),
                tags$p(info$Gene_Description))
      }
    )
    tier <- info$Tier
  }

  showModal(modalDialog(
    title = tagList(
      tags$span(symbol, style = "font-weight:700;font-size:1.2rem;"),
      if (!is.null(tier) && !is.na(tier))
        tags$span(tier, class = "badge bg-secondary",
                  style = "margin-left:8px;vertical-align:middle;")
    ),
    body,
    easyClose = TRUE,
    footer = modalButton("Close"),
    size = "l"
  ))
}

# Load + clean + annotate tier in one step.
load_annotated <- function(path) {
  annotate_tier(load_variants(path), TIER_DF)
}

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------
ui <- page_sidebar(
  title = "MacTel Variant Explorer",
  theme = bs_theme(version = 5, bootswatch = "flatly",
                   primary = "#1F4E79"),

  sidebar = sidebar(
    width = 320,
    title = "Filters",

    accordion(
      open = c("Data", "Core filters"),

      accordion_panel(
        "Data", icon = bsicons::bs_icon("database"),
        helpText(textOutput("data_source_label")),
        fileInput("upload", "Upload a Cavalier CSV",
                  accept = c(".csv"), buttonLabel = "Browse…"),
        actionButton("reset_filters", "Reset all filters",
                     class = "btn-outline-secondary btn-sm w-100")
      ),

      accordion_panel(
        "Core filters", icon = bsicons::bs_icon("funnel"),
        checkboxGroupInput("tier", "Tier",
                           choices = NULL, inline = TRUE),
        selectizeInput("genes", "Gene(s)", choices = NULL, multiple = TRUE,
                       options = list(placeholder = "All genes")),
        checkboxGroupInput("impact", "VEP impact",
                           choices = IMPACT_LEVELS, selected = IMPACT_LEVELS,
                           inline = TRUE),
        checkboxGroupInput("type", "Variant type",
                           choices = TYPE_LEVELS, selected = TYPE_LEVELS,
                           inline = TRUE)
      ),

      accordion_panel(
        "ClinVar & scores", icon = bsicons::bs_icon("sliders"),
        checkboxGroupInput("clnsig", "ClinVar class",
                           choices = CLNSIG_LEVELS, selected = CLNSIG_LEVELS),
        sliderInput("cadd", "CADD ≥",
                    min = 0, max = 60, value = 0, step = 1),
        sliderInput("revel", "REVEL ≥",
                    min = 0, max = 1, value = 0, step = 0.05),
        sliderInput("gnomad", "gnomAD AF ≤ (log10)",
                    min = -6, max = 0, value = 0, step = 0.5,
                    post = "")
      ),

      accordion_panel(
        "Priority flags", icon = bsicons::bs_icon("star"),
        helpText("Used by the 'Priority variants' tab only. A flag is set ",
                 "when a variant is ClinVar P/LP, HIGH impact, or CADD ≥ ",
                 "the threshold below. The tab keeps variants meeting at ",
                 "least the chosen number of flags."),
        sliderInput("min_flags", "Min. priority flags", min = 1, max = 3,
                    value = 1, step = 1),
        sliderInput("priority_cadd", "Flag: CADD ≥", min = 0, max = 60,
                    value = 20, step = 1)
      )
    )
  ),

  # value boxes (compact)
  layout_columns(
    fill = FALSE,
    value_box("Variants", textOutput("vb_variants"),
              showcase = bsicons::bs_icon("file-earmark-text"),
              theme = "primary", max_height = "92px"),
    value_box("Genes", textOutput("vb_genes"),
              showcase = bsicons::bs_icon("diagram-3"),
              theme = "secondary", max_height = "92px"),
    value_box("Samples", textOutput("vb_samples"),
              showcase = bsicons::bs_icon("people"),
              theme = "info", max_height = "92px"),
    value_box("ClinVar P/LP", textOutput("vb_plp"),
              showcase = bsicons::bs_icon("exclamation-triangle"),
              theme = "danger", max_height = "92px")
  ),

  navset_card_tab(
    id = "main_tabs",

    nav_panel(
      "Overview",
      layout_columns(
        col_widths = c(4, 4, 4),
        card(card_header("VEP impact"), plotOutput("p_impact", height = 280)),
        card(card_header("Variant type"), plotOutput("p_type", height = 280)),
        card(card_header("Inheritance"), plotOutput("p_inherit", height = 280))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("CADD distribution"), plotOutput("p_cadd", height = 320)),
        card(card_header("ClinVar classification"), plotOutput("p_clnsig", height = 320))
      ),
      card(card_header("Top genes by number of samples"),
           plotOutput("p_genes", height = 500))
    ),

    nav_panel(
      "Variant table",
      card(
        card_header(
          "Filtered variants",
          tags$span(bsicons::bs_icon("info-circle"),
                    " click a row for variant detail + protein lollipop",
                    class = "text-muted small ms-2"),
          downloadButton("dl_table", "Download CSV",
                         class = "btn-sm btn-primary float-end")
        ),
        DT::DTOutput("variant_table")
      )
    ),

    nav_panel(
      "Score scatter",
      card(
        card_header("CADD vs REVEL — hover for variant detail"),
        plotly::plotlyOutput("scatter", height = 600)
      )
    ),

    nav_panel(
      "Priority variants",
      card(
        card_header(
          "Priority variants (flag-filtered)",
          tags$span(bsicons::bs_icon("info-circle"),
                    " click a row for variant detail + protein lollipop",
                    class = "text-muted small ms-2"),
          downloadButton("dl_priority", "Download CSV",
                         class = "btn-sm btn-primary float-end")
        ),
        layout_columns(
          col_widths = c(7, 5),
          DT::DTOutput("priority_table"),
          plotOutput("p_priority_genes", height = 500)
        )
      )
    ),

    nav_panel(
      "Gene summary",
      card(
        card_header(
          "Per-gene summary",
          tags$span(bsicons::bs_icon("info-circle"),
                    " click a row for the gene description",
                    class = "text-muted small ms-2"),
          downloadButton("dl_genes", "Download CSV",
                         class = "btn-sm btn-primary float-end")
        ),
        DT::DTOutput("gene_table")
      )
    ),

    nav_panel(
      "Sample explorer",
      layout_sidebar(
        sidebar = sidebar(
          width = 280, position = "left",
          selectizeInput("sample_pick", "Select a sample",
                         choices = NULL, multiple = FALSE,
                         options = list(placeholder = "Start typing a sample ID…")),
          helpText("Shows every variant carried by the chosen sample, ",
                   "respecting the global filters in the left panel.")
        ),
        card(
          card_header(
            textOutput("sample_header"),
            tags$span(bsicons::bs_icon("info-circle"),
                      " click a row for variant detail + protein lollipop",
                      class = "text-muted small ms-2"),
            downloadButton("dl_sample", "Download CSV",
                           class = "btn-sm btn-primary float-end")
          ),
          DT::DTOutput("sample_table")
        )
      )
    )
  )
)

# -----------------------------------------------------------------------------
# SERVER
# -----------------------------------------------------------------------------
server <- function(input, output, session) {

  # ---- raw data (reactive on upload) ----------------------------------------
  raw <- reactiveVal(NULL)
  src_label <- reactiveVal("")

  observe({
    df <- load_annotated(startup_path)
    raw(df)
    src_label(sprintf("%s  (%d variants)",
                      basename(startup_path), nrow(df)))
  })

  observeEvent(input$upload, {
    req(input$upload)
    df <- tryCatch(load_annotated(input$upload$datapath),
                   error = function(e) {
                     showNotification(paste("Could not load file:", e$message),
                                      type = "error", duration = 8)
                     NULL
                   })
    if (!is.null(df)) {
      raw(df)
      src_label(sprintf("%s  (%d variants, uploaded)",
                        input$upload$name, nrow(df)))
    }
  })

  output$data_source_label <- renderText(paste("Source:", src_label()))

  # ---- populate dynamic filter choices when data changes --------------------
  observeEvent(raw(), {
    df <- raw()
    tiers <- sort(unique(df$Tier))
    updateCheckboxGroupInput(session, "tier",
                             choices = tiers, selected = tiers, inline = TRUE)
    updateSelectizeInput(session, "genes",
                         choices = sort(unique(df$SYMBOL)), server = TRUE)
    updateSelectizeInput(session, "sample_pick",
                         choices = sort(unique(df$family_id)), server = TRUE)
    mx <- ceiling(max(df$CADD, na.rm = TRUE))
    updateSliderInput(session, "cadd", max = mx, value = 0)
    updateSliderInput(session, "priority_cadd", max = mx)
  })

  observeEvent(input$reset_filters, {
    df <- raw(); req(df)
    updateCheckboxGroupInput(session, "tier",
                             selected = sort(unique(df$Tier)))
    updateSelectizeInput(session, "genes", selected = character(0))
    updateCheckboxGroupInput(session, "impact", selected = IMPACT_LEVELS)
    updateCheckboxGroupInput(session, "type", selected = TYPE_LEVELS)
    updateCheckboxGroupInput(session, "clnsig", selected = CLNSIG_LEVELS)
    updateSliderInput(session, "cadd", value = 0)
    updateSliderInput(session, "revel", value = 0)
    updateSliderInput(session, "gnomad", value = 0)
    updateSliderInput(session, "min_flags", value = 0)
  })

  # ---- core filtered dataset ------------------------------------------------
  filtered <- reactive({
    df <- raw(); req(df)

    if (length(input$tier) > 0)
      df <- dplyr::filter(df, Tier %in% input$tier)
    if (length(input$genes) > 0)
      df <- dplyr::filter(df, SYMBOL %in% input$genes)
    if (length(input$impact) > 0)
      df <- dplyr::filter(df, as.character(IMPACT) %in% input$impact)
    if (length(input$type) > 0)
      df <- dplyr::filter(df, as.character(TYPE) %in% input$type)
    if (length(input$clnsig) > 0)
      df <- dplyr::filter(df, as.character(CLNSIG_clean) %in% input$clnsig)

    df <- dplyr::filter(df, is.na(CADD) | CADD >= input$cadd)
    if (input$revel > 0)
      df <- dplyr::filter(df, is.na(REVEL) | REVEL >= input$revel)
    if (input$gnomad < 0) {
      thr <- 10^input$gnomad
      df <- dplyr::filter(df, is.na(gnomad_AF) | gnomad_AF <= thr)
    }

    # priority flags
    df <- df %>%
      dplyr::mutate(
        flag_clinvar = is_pathLP,
        flag_high    = IMPACT == "HIGH",
        flag_cadd    = !is.na(CADD) & CADD >= input$priority_cadd,
        n_flags      = as.integer(flag_clinvar) + as.integer(flag_high) +
                       as.integer(flag_cadd),
        why_prioritised = purrr_paste(flag_clinvar, flag_high, flag_cadd,
                                      CLNSIG_clean, CADD, IMPACT,
                                      input$priority_cadd)
      )
    df
  })

  # ---- priority subset (Priority variants tab) ------------------------------
  # Flags are computed in filtered(); here we keep only rows meeting at least
  # the chosen number of flags (min 1, so this is never the full table).
  priority <- reactive({
    filtered() %>%
      dplyr::filter(n_flags >= input$min_flags) %>%
      dplyr::arrange(dplyr::desc(n_flags), dplyr::desc(CADD))
  })

  # ---- value boxes ----------------------------------------------------------
  output$vb_variants <- renderText(format(nrow(filtered()), big.mark = ","))
  output$vb_genes    <- renderText(dplyr::n_distinct(filtered()$SYMBOL))
  output$vb_samples  <- renderText(dplyr::n_distinct(filtered()$family_id))
  output$vb_plp      <- renderText(sum(filtered()$is_pathLP))

  # ---- overview plots -------------------------------------------------------
  output$p_impact  <- renderPlot(plot_impact(filtered()))
  output$p_type    <- renderPlot(plot_type(filtered()))
  output$p_inherit <- renderPlot(plot_inheritance(filtered()))
  output$p_cadd    <- renderPlot(plot_cadd(filtered(), input$priority_cadd))
  output$p_clnsig  <- renderPlot(plot_clnsig(filtered()))
  output$p_genes   <- renderPlot(plot_top_genes(filtered(), 25))

  # ---- scatter --------------------------------------------------------------
  output$scatter <- plotly::renderPlotly({
    p <- plot_score_scatter(filtered())
    validate(need(!is.null(p),
                  "Need ≥3 variants with both CADD and REVEL scores."))
    plotly::ggplotly(p, tooltip = "text")
  })

  # ---- display-table builder ------------------------------------------------
  display_cols <- function(df) {
    df %>%
      dplyr::transmute(
        Gene = SYMBOL,
        Tier = Tier,
        Sample = family_id,
        Variant = paste0(CHROM, ":", POS, " ", REF, ">", ALT),
        HGVSc, HGVSp = HGVSp_short,
        Impact = IMPACT, Type = TYPE,
        CADD = round(CADD, 1),
        REVEL = round(REVEL, 3),
        AlphaMissense = am_class,
        SpliceAI = round(SpliceAI_max, 3),
        ClinVar = CLNSIG_clean,
        gnomAD_AF = signif(gnomad_AF, 3),
        Inheritance = inheritance
      )
  }

  output$variant_table <- DT::renderDT({
    DT::datatable(display_cols(filtered()),
                  filter = "top", rownames = FALSE,
                  selection = "single",
                  extensions = "Buttons",
                  options = list(pageLength = 25, scrollX = TRUE,
                                 dom = "Bfrtip", buttons = c("copy", "csv"))) %>%
      DT::formatStyle("ClinVar",
                      backgroundColor = DT::styleEqual(
                        c("Pathogenic", "Pathogenic/Likely_pathogenic",
                          "Likely_pathogenic"),
                        c("#FDDEDE", "#F7D6F7", "#FDE8D8"))) %>%
      DT::formatStyle("Impact",
                      backgroundColor = DT::styleEqual("HIGH", "#FDDEDE"))
  })

  # ---- priority table & plot ------------------------------------------------
  output$priority_table <- DT::renderDT({
    tbl <- priority() %>%
      dplyr::transmute(
        Gene = SYMBOL, Tier = Tier, Sample = family_id,
        HGVSc, HGVSp = HGVSp_short,
        Impact = IMPACT, Type = TYPE, CADD = round(CADD, 1),
        ClinVar = CLNSIG_clean, Flags = n_flags,
        `Why prioritised` = why_prioritised)
    DT::datatable(tbl, filter = "top", rownames = FALSE,
                  selection = "single",
                  options = list(pageLength = 15, scrollX = TRUE)) %>%
      DT::formatStyle("Flags", fontWeight = "bold",
                      background = DT::styleColorBar(c(0, 3), "#9ec5fe")) %>%
      DT::formatStyle("Impact",
                      backgroundColor = DT::styleEqual("HIGH", "#FDDEDE"))
  })

  output$p_priority_genes <- renderPlot({
    d <- priority()
    validate(need(nrow(d) > 0, "No variants meet the chosen number of flags."))
    plot_top_genes(d, 20)
  })

  # ---- gene summary ---------------------------------------------------------
  gene_summary <- reactive({
    filtered() %>%
      dplyr::group_by(SYMBOL, Tier) %>%
      dplyr::summarise(
        Variants    = dplyr::n(),
        Samples     = dplyr::n_distinct(family_id),
        `P/LP`      = sum(is_pathLP),
        `HIGH`      = sum(IMPACT == "HIGH"),
        CADD_max    = round(max(CADD, na.rm = TRUE), 1),
        REVEL_max   = round(suppressWarnings(max(REVEL, na.rm = TRUE)), 3),
        Types       = paste(sort(unique(as.character(TYPE))), collapse = "/"),
        .groups = "drop") %>%
      dplyr::mutate(REVEL_max = ifelse(is.infinite(REVEL_max), NA, REVEL_max)) %>%
      dplyr::arrange(dplyr::desc(Samples), dplyr::desc(CADD_max))
  })

  output$gene_table <- DT::renderDT({
    DT::datatable(gene_summary(), rownames = FALSE, filter = "top",
                  selection = "single",
                  options = list(pageLength = 25, scrollX = TRUE))
  })

  # ---- sample explorer ------------------------------------------------------
  sample_data <- reactive({
    req(input$sample_pick)
    filtered() %>%
      dplyr::filter(family_id == input$sample_pick) %>%
      dplyr::arrange(dplyr::desc(is_pathLP), dplyr::desc(CADD))
  })

  output$sample_header <- renderText({
    if (is.null(input$sample_pick) || input$sample_pick == "")
      "Select a sample to see its variants"
    else sprintf("Variants in sample %s", input$sample_pick)
  })

  output$sample_table <- DT::renderDT({
    validate(need(!is.null(input$sample_pick) && input$sample_pick != "",
                  "Choose a sample from the dropdown."))
    d <- sample_data()
    validate(need(nrow(d) > 0,
                  "This sample has no variants under the current filters."))
    tbl <- d %>%
      dplyr::transmute(
        Gene = SYMBOL, Tier = Tier,
        Variant = paste0(CHROM, ":", POS, " ", REF, ">", ALT),
        HGVSc, HGVSp = HGVSp_short,
        Impact = IMPACT, Type = TYPE, CADD = round(CADD, 1),
        REVEL = round(REVEL, 3), AlphaMissense = am_class,
        ClinVar = CLNSIG_clean, gnomAD_AF = signif(gnomad_AF, 3),
        Inheritance = inheritance)
    DT::datatable(tbl, rownames = FALSE,
                  selection = "single",
                  options = list(pageLength = 25, scrollX = TRUE)) %>%
      DT::formatStyle("ClinVar",
                      backgroundColor = DT::styleEqual(
                        c("Pathogenic", "Pathogenic/Likely_pathogenic",
                          "Likely_pathogenic"),
                        c("#FDDEDE", "#F7D6F7", "#FDE8D8"))) %>%
      DT::formatStyle("Impact",
                      backgroundColor = DT::styleEqual("HIGH", "#FDDEDE"))
  })

  # ---- click a table row -> gene description modal --------------------------
  # Row-selection indices from DT refer to the data in its original order,
  # so we re-derive each table's data frame to map row -> gene symbol.
  observeEvent(input$variant_table_rows_selected, {
    sel <- input$variant_table_rows_selected
    df  <- display_cols(filtered())
    if (length(sel) && sel <= nrow(df)) show_gene_modal(df$Gene[sel])
    DT::dataTableProxy("variant_table") %>% DT::selectRows(NULL)
  })

  observeEvent(input$priority_table_rows_selected, {
    sel <- input$priority_table_rows_selected
    df  <- priority()
    if (length(sel) && sel <= nrow(df)) show_gene_modal(df$SYMBOL[sel])
    DT::dataTableProxy("priority_table") %>% DT::selectRows(NULL)
  })

  observeEvent(input$gene_table_rows_selected, {
    sel <- input$gene_table_rows_selected
    df  <- gene_summary()
    if (length(sel) && sel <= nrow(df)) show_gene_modal(df$SYMBOL[sel])
    DT::dataTableProxy("gene_table") %>% DT::selectRows(NULL)
  })

  observeEvent(input$sample_table_rows_selected, {
    sel <- input$sample_table_rows_selected
    df  <- sample_data()
    if (length(sel) && sel <= nrow(df)) show_gene_modal(df$SYMBOL[sel])
    DT::dataTableProxy("sample_table") %>% DT::selectRows(NULL)
  })

  # ---- downloads ------------------------------------------------------------
  output$dl_table <- downloadHandler(
    filename = function() sprintf("filtered_variants_%s.csv", Sys.Date()),
    content  = function(file) readr::write_csv(display_cols(filtered()), file)
  )
  output$dl_priority <- downloadHandler(
    filename = function() sprintf("priority_variants_%s.csv", Sys.Date()),
    content  = function(file) {
      d <- priority() %>%
        dplyr::select(SYMBOL, Tier, family_id, CHROM, POS, REF, ALT,
                      HGVSc, HGVSp_short, IMPACT, TYPE, CADD, REVEL,
                      am_class, SpliceAI_max, CLNSIG_clean, gnomad_AF,
                      inheritance, flag_clinvar, flag_high, flag_cadd,
                      n_flags, why_prioritised) %>%
        dplyr::rename(Sample = family_id)
      readr::write_csv(d, file)
    }
  )
  output$dl_genes <- downloadHandler(
    filename = function() sprintf("gene_summary_%s.csv", Sys.Date()),
    content  = function(file) readr::write_csv(gene_summary(), file)
  )
  output$dl_sample <- downloadHandler(
    filename = function() sprintf("sample_%s_variants_%s.csv",
                                  input$sample_pick %||% "none", Sys.Date()),
    content  = function(file) {
      d <- sample_data() %>%
        dplyr::select(family_id, SYMBOL, Tier, CHROM, POS, REF, ALT,
                      HGVSc, HGVSp_short, IMPACT, TYPE, CADD, REVEL,
                      am_class, SpliceAI_max, CLNSIG_clean, gnomad_AF,
                      inheritance) %>%
        dplyr::rename(Sample = family_id)
      readr::write_csv(d, file)
    }
  )
}

# helper used inside the reactive (vectorised "why prioritised" text) ---------
purrr_paste <- function(fc, fh, fca, clnsig, cadd, impact, cadd_thr) {
  vapply(seq_along(fc), function(i) {
    parts <- character(0)
    if (isTRUE(fc[i]))  parts <- c(parts, gsub("_", " ", as.character(clnsig[i])))
    if (isTRUE(fh[i]))  parts <- c(parts, "HIGH impact")
    if (isTRUE(fca[i])) parts <- c(parts, paste0("CADD ", round(cadd[i], 1)))
    paste(parts, collapse = " | ")
  }, character(1))
}

shinyApp(ui, server)

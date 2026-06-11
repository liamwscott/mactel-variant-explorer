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

DEFAULT_REAL    <- file.path(app_dir, "data", "candidate_variants_tier1_2.csv")
DEFAULT_EXAMPLE <- file.path(app_dir, "data", "example_variants.csv")
startup_path    <- if (file.exists(DEFAULT_REAL)) DEFAULT_REAL else DEFAULT_EXAMPLE

# Gene info (Tier + descriptions; gene symbols only, no participant data) ------
GENE_INFO_PATH <- file.path(app_dir, "data", "gene_info.tsv")
GENE_INFO      <- load_gene_info(GENE_INFO_PATH)

# Pfam protein domains, pre-fetched offline (gene symbols only, no patient data).
PROTEIN_DOMAINS <- load_protein_domains(file.path(app_dir, "data",
                                                  "protein_domains.tsv"))

# Per-sample info (case/control status + data-group flags). Prefer the real
# sheet when present (git-ignored); otherwise the de-identified example.
SAMPLE_INFO_REAL    <- file.path(app_dir, "data", "all_samples_fixed.txt")
SAMPLE_INFO_EXAMPLE <- file.path(app_dir, "data", "example_sample_info.tsv")
SAMPLE_INFO <- load_sample_info(
  if (file.exists(SAMPLE_INFO_REAL)) SAMPLE_INFO_REAL else SAMPLE_INFO_EXAMPLE)

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
    footer = tagList(
      tags$button(
        tagList(bsicons::bs_icon("table"), " List all variants"),
        class = "btn btn-primary",
        onclick = sprintf(
          "Shiny.setInputValue('gene_view_variants','%s',{priority:'event'});return false;",
          gsub("'", "\\\\'", symbol)
        )
      ),
      tags$button(
        tagList(bsicons::bs_icon("graph-up"), " View variants on protein"),
        class = "btn btn-outline-primary",
        onclick = sprintf(
          "Shiny.setInputValue('gene_view_lollipop','%s',{priority:'event'});return false;",
          gsub("'", "\\\\'", symbol)
        )
      ),
      modalButton("Close")
    ),
    size = "l"
  ))
}

# Load + clean + annotate tier in one step.
load_annotated <- function(path) {
  annotate_tier(load_variants(path), TIER_DF)
}

# Short descriptions shown as hover tooltips on table column headers. Keyed by
# the displayed column label; covers every column used across the app's tables,
# so the same callback can be reused everywhere (unmatched headers are ignored).
COLUMN_TIPS <- c(
  "Gene"          = "Gene symbol — click to open its description",
  "SYMBOL"        = "Gene symbol — click to open its description",
  "Tier"          = "Curated MacTel gene tier (Tier 1 = strongest evidence)",
  "Sample"        = "Sample / individual ID — click to open in the sample explorer",
  "Variant"       = "Genomic change CHROM:POS REF>ALT — click for the protein lollipop",
  "HGVSc"         = "Coding-DNA change in HGVS c. notation",
  "HGVSp"         = "Protein change in HGVS p. notation",
  "Impact"        = "VEP-predicted consequence severity (HIGH / MODERATE / LOW / MODIFIER)",
  "Type"          = "Variant class (LOF / SPLICING / MISSENSE / OTHER)",
  "CADD"          = "CADD deleteriousness score (PHRED-scaled; >20 ~ top 1% most deleterious)",
  "REVEL"         = "REVEL missense pathogenicity score (0-1; higher = more damaging)",
  "AlphaMissense" = "AlphaMissense class (likely benign / ambiguous / likely pathogenic)",
  "SpliceAI"      = "Maximum SpliceAI delta score (0-1; higher = stronger predicted splice effect)",
  "ClinVar"       = "ClinVar clinical-significance classification",
  "gnomAD_AF"     = "gnomAD population allele frequency",
  "Inheritance"   = "Inheritance pattern observed for this variant",
  "Flags"         = "Number of priority flags met (ClinVar P/LP, HIGH impact, high CADD)",
  "Why prioritised" = "Which priority criteria this variant meets",
  "Variants"      = "Number of variants in this gene under the current filters",
  "Samples"       = "Number of distinct samples carrying a variant in this gene",
  "P/LP"          = "Count of Pathogenic / Likely-pathogenic variants (ClinVar)",
  "HIGH"          = "Count of HIGH-impact variants (VEP)",
  "CADD_max"      = "Highest CADD score among this gene's variants",
  "REVEL_max"     = "Highest REVEL score among this gene's variants",
  "Types"         = "Variant classes present in this gene"
)

# headerCallback that attaches the COLUMN_TIPS as native title= tooltips. Matches
# by header text, so it coexists with the filter row and works on any table.
header_tips_cb <- function() {
  tips_json <- jsonlite::toJSON(as.list(COLUMN_TIPS), auto_unbox = TRUE)
  DT::JS(sprintf(
    "function(thead, data, start, end, display) {
       var tips = %s;
       $(thead).find('th').each(function() {
         var t = $(this).text().trim();
         if (tips[t]) { $(this).attr('title', tips[t]); }
       });
     }", tips_json))
}

# -----------------------------------------------------------------------------
# UI
# -----------------------------------------------------------------------------
# App theme. Pagination is forced to the sidebar blue (flatly otherwise renders
# the page-number controls in its default green) by compiling the override into
# the theme stylesheet, which is more reliable than injecting <style> tags.
app_theme <- bs_theme(version = 5, bootswatch = "flatly", primary = "#1F4E79") |>
  bslib::bs_add_variables(
    "pagination-color"                 = "#1F4E79",
    "pagination-hover-color"           = "#ffffff",
    "pagination-hover-bg"              = "#1F4E79",
    "pagination-hover-border-color"    = "#1F4E79",
    "pagination-active-bg"             = "#1F4E79",
    "pagination-active-border-color"   = "#1F4E79",
    .where = "declarations"
  ) |>
  bslib::bs_add_rules("
    .dataTables_wrapper .dataTables_paginate .paginate_button { color: #1F4E79 !important; }
    .dataTables_wrapper .dataTables_paginate .paginate_button.current,
    .dataTables_wrapper .dataTables_paginate .paginate_button.current:hover,
    .dataTables_wrapper .dataTables_paginate .paginate_button:hover {
      color: #ffffff !important;
      background: #1F4E79 !important;
      background-image: none !important;
      border: 1px solid #1F4E79 !important;
    }
    .page-link { color: #1F4E79 !important; }
    .page-item.active .page-link {
      background-color: #1F4E79 !important;
      border-color: #1F4E79 !important;
      color: #ffffff !important;
    }

    /* Tighter tab bar — smaller text/padding so all tabs fit on one row. */
    .nav-tabs .nav-link {
      padding: 0.3rem 0.65rem !important;
      font-size: 0.82rem !important;
    }

    /* Smaller value-box icons and text in the header. */
    .bslib-value-box .value-box-title { font-size: 0.78rem !important; margin-bottom: 0 !important; }
    .bslib-value-box .value-box-value { font-size: 1.05rem !important; }
    .bslib-value-box .value-box-showcase { padding: 0.25rem 0.5rem !important; }
    .bslib-value-box .value-box-showcase svg,
    .bslib-value-box .value-box-showcase .bi {
      width: 1.4rem !important;
      height: 1.4rem !important;
      font-size: 1.4rem !important;
    }
  ")

# Glossary helper for the "Start here" tab. Renders one collapsible entry with
# a bold one-line summary followed by a plain-English explanation.
gloss <- function(title, summary, ...) {
  accordion_panel(
    title,
    tags$p(class = "mb-1", tags$strong(summary)),
    tags$div(class = "text-body-secondary small", ...)
  )
}

# Reusable DNA-helix logo (used in the title bar and the Genes value box).
dna_icon <- function(size = 24) HTML(paste0(
  '<svg xmlns="http://www.w3.org/2000/svg" width="', size, '" height="', size, '" ',
  'viewBox="0 0 24 24" fill="none" stroke="currentColor" stroke-width="2" ',
  'stroke-linecap="round" stroke-linejoin="round" style="vertical-align:middle;">',
  '<path d="m10 16 1.5 1.5"/><path d="m14 8-1.5-1.5"/>',
  '<path d="M15 2c-1.798 1.998-2.518 3.995-2.807 5.993"/>',
  '<path d="m16.5 10.5 1 1"/><path d="m17 6-2.891-2.891"/>',
  '<path d="M2 15c6.667-6 13.333 0 20-6"/><path d="m20 9 .891.891"/>',
  '<path d="M3.109 14.109 4 15"/><path d="m6.5 12.5 1 1"/>',
  '<path d="m7 18 2.891 2.891"/>',
  '<path d="M9 22c1.798-1.998 2.518-3.995 2.807-5.993"/></svg>'))

# A numbered "quick start" step for the landing page.
landing_step <- function(n, title, body) div(
  class = "d-flex mb-3",
  div(class = "flex-shrink-0 d-flex align-items-center justify-content-center",
      style = paste("width:30px;height:30px;border-radius:50%;background:#1F4E79;",
                    "color:#fff;font-weight:600;font-size:0.9rem;"),
      n),
  div(class = "ms-3",
      tags$div(tags$strong(title)),
      tags$div(class = "small text-body-secondary", body))
)

# A tab entry for the landing-page "tabs" list, with its matching icon.
tab_item <- function(icon, name, desc) tags$li(
  class = "mb-2",
  tags$span(class = "text-primary", bsicons::bs_icon(icon)),
  tags$strong(paste0(" ", name)), desc)

ui <- function(request) page_sidebar(
  title = tags$span(
    class = "d-inline-flex align-items-center",
    dna_icon(32),
    tags$span("MacTel Variant Explorer", class = "ms-2",
              style = "font-size:1.75rem;font-weight:400;letter-spacing:0.2px;")
  ),
  theme = app_theme,
  # Non-fillable so cards keep their real heights and the page scrolls,
  # rather than squeezing every plot into a single viewport.
  fillable = FALSE,

  sidebar = sidebar(
    width = 320,
    title = "Filters",

    actionButton("reset_filters",
                 tagList(bsicons::bs_icon("arrow-counterclockwise"),
                         " Reset all filters"),
                 class = "btn-outline-primary btn-sm w-100 mb-2"),

    accordion(
      open = c("Data", "Core filters"),

      accordion_panel(
        "Data", icon = bsicons::bs_icon("database"),
        helpText(textOutput("data_source_label")),
        fileInput("upload", "Upload a Cavalier CSV",
                  accept = c(".csv"), buttonLabel = "Browse…")
      ),

      accordion_panel(
        "Core filters", icon = bsicons::bs_icon("funnel"),
        checkboxGroupInput("sample_group", "Sample group",
                           choices = c("MacTel", "HSAN1", "Controls"),
                           selected = "MacTel", inline = TRUE),
        checkboxGroupInput("tier", "Gene Tier",
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
    ),

    actionButton("share_link",
                 tagList(bsicons::bs_icon("link-45deg"), " Copy share link"),
                 class = "btn-outline-primary btn-sm w-100 mt-3",
                 title = paste("Get a link that reproduces the current filters",
                               "and tab — paste it to a colleague or bookmark it."))
  ),

  # value boxes (compact)
  layout_columns(
    fill = FALSE,
    value_box("Variants", textOutput("vb_variants"),
              showcase = bsicons::bs_icon("file-earmark-text"),
              showcase_layout = bslib::showcase_left_center(width = "3rem"),
              theme = "primary", max_height = "80px"),
    value_box("Genes", textOutput("vb_genes"),
              showcase = dna_icon(),
              showcase_layout = bslib::showcase_left_center(width = "3rem"),
              theme = "secondary", max_height = "80px"),
    value_box("Samples", textOutput("vb_samples"),
              showcase = bsicons::bs_icon("people"),
              showcase_layout = bslib::showcase_left_center(width = "3rem"),
              theme = "info", max_height = "80px"),
    value_box("ClinVar P/LP", textOutput("vb_plp"),
              showcase = bsicons::bs_icon("exclamation-triangle"),
              showcase_layout = bslib::showcase_left_center(width = "3rem"),
              theme = "danger", max_height = "80px")
  ),

  navset_card_tab(
    id = "main_tabs",

    nav_panel(
      "Start here",
      icon = bsicons::bs_icon("compass"),
      # Hero banner
      div(
        class = "p-4 mb-3 rounded-3 shadow-sm",
        style = "background:linear-gradient(135deg,#1F4E79 0%,#3A7CA5 100%);color:#fff;",
        div(
          class = "d-flex align-items-center",
          tags$span(class = "d-inline-flex", style = "opacity:.95;", dna_icon(54)),
          div(
            class = "ms-3",
            tags$h2("Welcome to the MacTel Variant Explorer",
                    class = "mb-1", style = "font-weight:600;font-size:1.6rem;"),
            tags$p(class = "mb-0", style = "opacity:.92;font-size:1.02rem;",
              "Explore, filter, and prioritise rare genetic variants from the ",
              "MacTel study — no command line required.")
          )
        )
      ),
      layout_columns(
        col_widths = c(7, 5),
        card(
          card_header(bsicons::bs_icon("info-circle"),
                      tags$strong(" How to use this app")),
          tags$p(class = "text-body-secondary small mb-3",
            "Each row in the data is a ", tags$strong("variant"),
            " (a single change in a person's DNA) seen in a ",
            tags$strong("sample"), " (one participant), in a MacTel gene."),
          landing_step(1, "Filter on the left",
            "Narrow by gene, sample group, predicted severity, or how rare the variant is. The counters up top update live."),
          landing_step(2, "Browse the tabs",
            "Each tab shows the filtered variants a different way — charts, tables, and plots."),
          landing_step(3, "Click anything blue",
            "Genes, variants, and samples are links that open detail views."),
          div(class = "alert alert-primary d-flex align-items-center mb-0 py-2",
              role = "alert",
              bsicons::bs_icon("hand-index-thumb"),
              tags$span(class = "ms-2 small",
                tags$strong("Tip: "), "click a ", tags$strong("gene"),
                " for its description, a ", tags$strong("variant"),
                " to see it on the protein, or a ", tags$strong("sample"),
                " to open that participant's profile."))
        ),
        card(
          card_header(bsicons::bs_icon("signpost-2"), tags$strong(" The tabs")),
          tags$ul(class = "list-unstyled mb-0 small",
            tab_item("bar-chart-line", "Overview",
                     " — summary charts of the variants currently filtered in."),
            tab_item("table", "Variant table",
                     " — every filtered variant in a searchable, sortable table."),
            tab_item("graph-up", "Score scatter",
                     " — CADD vs REVEL, to spot variants high on both."),
            tab_item("star-fill", "Priority variants",
                     " — the strongest candidates, with a plain-English reason."),
            tab_item("card-list", "Gene summary",
                     " — one row per gene, rolling up its variants."),
            tab_item("person-lines-fill", "Sample explorer",
                     " — everything for a single participant.")
          )
        )
      ),
      card(
        card_header(bsicons::bs_icon("book"),
                    tags$strong(" Glossary"),
                    tags$span(" — what do these terms and scores mean?",
                              class = "text-body-secondary")),
        tags$p(class = "text-body-secondary small",
          "Click any term to expand it. The numeric cut-offs below are common ",
          "rules of thumb, not hard rules — always interpret a variant in ",
          "context."),
        accordion(
          open = FALSE,
          gloss("CADD", "How damaging a variant is predicted to be (any variant type).",
                "Scaled 0–99. Higher means more likely to be harmful. As a guide, ",
                "a score of ", tags$strong("20"), " puts a variant in the top 1% ",
                "most deleterious in the genome, and ", tags$strong("30"),
                " in the top 0.1%."),
          gloss("REVEL", "Likelihood that a missense change is disease-causing.",
                "A score from 0 to 1 for ", tags$strong("missense"),
                " variants (one amino acid swapped for another). Higher means ",
                "more likely pathogenic; values above ~0.5 are suggestive and ",
                "above ~0.75 are stronger evidence."),
          gloss("AlphaMissense", "Google DeepMind's AI prediction for missense changes.",
                "Gives each missense variant a score (0–1) and a class: ",
                tags$strong("likely_benign"), ", ", tags$strong("ambiguous"),
                ", or ", tags$strong("likely_pathogenic"),
                ". Likely-pathogenic calls are highlighted in the tables."),
          gloss("SpliceAI", "Whether a variant is predicted to disrupt splicing.",
                "Splicing is how the cell stitches a gene's coding pieces ",
                "together; disrupting it can break the protein. Scored 0–1: ",
                "≥0.2 possible, ≥0.5 likely, ≥0.8 high-confidence splice effect."),
          gloss("gnomAD allele frequency (AF)", "How common the variant is in the general population.",
                "From the gnomAD reference database of >100,000 people. ",
                tags$strong("Lower is rarer."), " Disease-causing variants for a ",
                "rare condition are usually very rare (e.g. AF below 0.001). The ",
                "slider uses log10, so −3 means AF ≤ 0.001."),
          gloss("ClinVar significance", "What clinical databases say about the variant.",
                "ClinVar is a public archive of variant interpretations: ",
                tags$strong("Pathogenic"), ", ", tags$strong("Likely pathogenic"),
                ", ", tags$strong("Uncertain significance (VUS)"), ", ",
                tags$strong("Conflicting"), ", ",
                tags$strong("Benign/Likely benign"), ", or ",
                tags$strong("Not in ClinVar"), " if it has never been submitted."),
          gloss("VEP impact", "A severity category for the variant's effect on the gene.",
                tags$strong("HIGH"), " (e.g. a premature stop, frameshift, or ",
                "splice-site change — likely to break the protein), ",
                tags$strong("MODERATE"), " (e.g. a missense change), ",
                tags$strong("LOW"), " (e.g. a silent change), and ",
                tags$strong("MODIFIER"), " (non-coding / regulatory regions)."),
          gloss("Variant type", "A simplified grouping of the change.",
                tags$strong("LOF"), " (loss of function — disables the gene), ",
                tags$strong("SPLICING"), ", ", tags$strong("MISSENSE"),
                " (amino-acid change), or ", tags$strong("OTHER"), "."),
          gloss("Inheritance", "The predicted way the variant was inherited.",
                "Based on the genotype pattern — for example de novo (new in the ",
                "child), dominant, recessive (two copies), compound heterozygous, ",
                "or X-linked. Shown as \"unknown\" when it can't be determined."),
          gloss("Gene tier", "How strong the evidence is that the gene is involved in MacTel.",
                tags$strong("Tier 1"), " genes have established, high-confidence ",
                "links to MacTel; ", tags$strong("Tier 2"),
                " are candidate genes with supporting but less definitive ",
                "evidence. ", tags$strong("Unassigned"),
                " means the gene isn't on the curated list."),
          gloss("Sample groups", "How participants are categorised.",
                tags$strong("MacTel"), " — diagnosed cases; ",
                tags$strong("HSAN1"), " — carries an HSAN1-causing variant ",
                "(in SPTLC1/SPTLC2); ", tags$strong("Controls"),
                " — participants who are neither."),
          gloss("HGVSc / HGVSp", "The standard names for a variant.",
                tags$strong("HGVSc"), " describes the change at the DNA/coding ",
                "level (starts with \"c.\"); ", tags$strong("HGVSp"),
                " describes the resulting protein change (starts with \"p.\").")
        )
      )
    ),

    nav_panel(
      "Overview",
      icon = bsicons::bs_icon("bar-chart-line"),
      layout_columns(
        col_widths = c(4, 4, 4),
        card(card_header("VEP impact"), plotOutput("p_impact", height = 340)),
        card(card_header("Variant type"), plotOutput("p_type", height = 340)),
        card(card_header("Inheritance"), plotOutput("p_inherit", height = 340))
      ),
      layout_columns(
        col_widths = c(6, 6),
        card(card_header("CADD distribution"), plotOutput("p_cadd", height = 380)),
        card(card_header("ClinVar classification"), plotOutput("p_clnsig", height = 380))
      ),
      card(card_header("Top genes by number of samples"),
           plotOutput("p_genes", height = 500))
    ),

    nav_panel(
      "Variant table",
      icon = bsicons::bs_icon("table"),
      card(
        card_header(
          "Filtered variants",
          tags$span(bsicons::bs_icon("info-circle"),
                    " click Gene for its description, Variant for the lollipop, ",
                    "or Sample to open the sample explorer",
                    class = "text-muted small ms-2"),
          downloadButton("dl_table", "Download CSV",
                         class = "btn-sm btn-primary float-end")
        ),
        DT::DTOutput("variant_table")
      )
    ),

    nav_panel(
      "Score scatter",
      icon = bsicons::bs_icon("graph-up"),
      card(
        card_header("CADD vs REVEL — hover for variant detail"),
        plotly::plotlyOutput("scatter", height = 600)
      )
    ),

    nav_panel(
      "Priority variants",
      icon = bsicons::bs_icon("star-fill"),
      card(
        card_header(
          "Priority variants (flag-filtered)",
          tags$span(bsicons::bs_icon("info-circle"),
                    " click Gene for its description, Variant for the lollipop, ",
                    "or Sample to open the sample explorer",
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
      icon = bsicons::bs_icon("card-list"),
      card(
        card_header(
          "Per-gene summary",
          tags$span(bsicons::bs_icon("info-circle"),
                    " click a gene for its description",
                    class = "text-muted small ms-2"),
          downloadButton("dl_genes", "Download CSV",
                         class = "btn-sm btn-primary float-end")
        ),
        DT::DTOutput("gene_table")
      )
    ),

    nav_panel(
      "Sample explorer",
      icon = bsicons::bs_icon("person-lines-fill"),
      layout_sidebar(
        sidebar = sidebar(
          width = 280, position = "left",
          selectizeInput("sample_pick", "Select a sample",
                         choices = NULL, multiple = FALSE,
                         options = list(placeholder = "Start typing a sample ID…")),
          helpText("Shows every variant carried by the chosen sample. ",
                   "The lower table ignores the global filters; the upper one ",
                   "respects them.")
        ),
        div(
          uiOutput("sample_tags"),
          card(
            card_header(
              "Filtered / prioritised variants",
              tags$span(bsicons::bs_icon("info-circle"),
                        " variants for this sample that pass the global filters",
                        class = "text-muted small ms-2")
            ),
            DT::DTOutput("sample_table_priority")
          ),
          card(
            card_header(
              textOutput("sample_header"),
              tags$span(bsicons::bs_icon("info-circle"),
                        " every variant this sample carries, ignoring filters",
                        class = "text-muted small ms-2"),
              downloadButton("dl_sample", "Download CSV",
                             class = "btn-sm btn-primary float-end")
            ),
            DT::DTOutput("sample_table_all")
          )
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
    updateCheckboxGroupInput(session, "sample_group", selected = "MacTel")
  })

  # ---- core filtered dataset ------------------------------------------------
  # All filters EXCEPT the sample-group filter. The sample explorer uses this so
  # an individual can be searched even if their diagnosis group is unticked.
  filtered_pre_group <- reactive({
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
      df <- dplyr::filter(df, !is.na(REVEL) & REVEL >= input$revel)
    if (input$gnomad < 0) {
      thr <- 10^input$gnomad
      df <- dplyr::filter(df, is.na(gnomad_AF) | gnomad_AF <= thr)
    }
    df
  })

  filtered <- reactive({
    df <- filtered_pre_group()

    # Sample-group filter: union of ticked groups (MacTel and HSAN1 may
    # overlap). No ticks shows nothing.
    if (!is.null(SAMPLE_INFO)) {
      sel <- input$sample_group
      allowed <- SAMPLE_INFO$family_id[
        (("MacTel"   %in% sel) & SAMPLE_INFO$is_mactel)  |
        (("HSAN1"    %in% sel) & SAMPLE_INFO$is_hsan1)   |
        (("Controls" %in% sel) & SAMPLE_INFO$is_control)]
      df <- dplyr::filter(df, family_id %in% allowed)
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

  # ---- clickable-cell link builders -----------------------------------------
  # Each renders an <a> that fires a Shiny input carrying the row's identity, so
  # Gene/Variant/Sample cells become independent click targets in the same row.
  # The JS-string payload is sorted-on by DT (prefix is constant), so columns
  # still sort by the embedded id.
  .jsesc <- function(x) gsub("'", "\\\\'", as.character(x))
  link_gene <- function(sym) {
    ifelse(is.na(sym) | sym == "", as.character(sym), sprintf(
      "<a href='#' onclick=\"Shiny.setInputValue('cell_gene','%s',{priority:'event'});return false;\">%s</a>",
      .jsesc(sym), sym))
  }
  link_sample <- function(fid) sprintf(
    "<a href='#' onclick=\"Shiny.setInputValue('cell_sample','%s',{priority:'event'});return false;\">%s</a>",
    .jsesc(fid), fid)
  link_variant <- function(chrom, pos, ref, alt) {
    label <- sprintf("%s:%s %s>%s", chrom, pos, ref, alt)
    key   <- sprintf("%s||%s||%s||%s", chrom, pos, ref, alt)
    sprintf(
      "<a href='#' onclick=\"Shiny.setInputValue('cell_variant','%s',{priority:'event'});return false;\">%s</a>",
      .jsesc(key), label)
  }

  # ---- display-table builder ------------------------------------------------
  display_cols <- function(df, links = TRUE) {
    if (links) {
      gene_col    <- link_gene(df$SYMBOL)
      sample_col  <- link_sample(df$family_id)
      variant_col <- link_variant(df$CHROM, df$POS, df$REF, df$ALT)
    } else {
      gene_col    <- df$SYMBOL
      sample_col  <- df$family_id
      variant_col <- sprintf("%s:%s %s>%s", df$CHROM, df$POS, df$REF, df$ALT)
    }
    df %>%
      dplyr::transmute(
        Gene = gene_col,
        Tier = Tier,
        Sample = sample_col,
        Variant = variant_col,
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
                  selection = "none", escape = FALSE,
                  extensions = "Buttons",
                  options = list(pageLength = 25, scrollX = TRUE,
                                 dom = "Bfrtip", buttons = c("copy", "csv"),
                                 headerCallback = header_tips_cb())) %>%
      DT::formatStyle("ClinVar",
                      backgroundColor = DT::styleEqual(
                        c("Pathogenic", "Pathogenic/Likely_pathogenic",
                          "Likely_pathogenic"),
                        c("#FDDEDE", "#F7D6F7", "#FDE8D8"))) %>%
      DT::formatStyle("AlphaMissense",
                      backgroundColor = DT::styleEqual(
                        "likely_pathogenic", "#FDE8D8")) %>%
      DT::formatStyle("Impact",
                      backgroundColor = DT::styleEqual("HIGH", "#FDDEDE"))
  })

  # ---- priority table & plot ------------------------------------------------
  output$priority_table <- DT::renderDT({
    tbl <- priority() %>%
      dplyr::transmute(
        Gene = link_gene(SYMBOL), Tier = Tier, Sample = link_sample(family_id),
        Variant = link_variant(CHROM, POS, REF, ALT),
        HGVSc, HGVSp = HGVSp_short,
        Impact = IMPACT, Type = TYPE, CADD = round(CADD, 1),
        ClinVar = CLNSIG_clean, Flags = n_flags,
        `Why prioritised` = why_prioritised)
    DT::datatable(tbl, filter = "top", rownames = FALSE,
                  selection = "none", escape = FALSE,
                  options = list(pageLength = 15, scrollX = TRUE,
                                 headerCallback = header_tips_cb())) %>%
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
    g <- gene_summary() %>% dplyr::mutate(SYMBOL = link_gene(SYMBOL))
    DT::datatable(g, rownames = FALSE, filter = "top",
                  selection = "none", escape = FALSE,
                  options = list(pageLength = 25, scrollX = TRUE,
                                 headerCallback = header_tips_cb()))
  })

  # ---- sample explorer ------------------------------------------------------
  # Prioritised view respects the global filters (but not the sample-group one);
  # the "all variants" view ignores every filter so the full carrier set shows.
  sample_data_priority <- reactive({
    req(input$sample_pick)
    filtered_pre_group() %>%
      dplyr::filter(family_id == input$sample_pick) %>%
      dplyr::arrange(dplyr::desc(is_pathLP), dplyr::desc(CADD))
  })

  sample_data_all <- reactive({
    req(input$sample_pick)
    raw() %>%
      dplyr::filter(family_id == input$sample_pick) %>%
      dplyr::arrange(dplyr::desc(is_pathLP), dplyr::desc(CADD))
  })

  output$sample_header <- renderText({
    if (is.null(input$sample_pick) || input$sample_pick == "")
      "Select a sample to see its variants"
    else sprintf("Variants in sample %s", input$sample_pick)
  })

  # Data-group tags for the picked sample (Mito_haplo deliberately excluded).
  output$sample_tags <- renderUI({
    if (is.null(SAMPLE_INFO) ||
        is.null(input$sample_pick) || input$sample_pick == "")
      return(NULL)
    row <- SAMPLE_INFO[SAMPLE_INFO$family_id == input$sample_pick, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)
    row <- row[1, ]

    # Diagnosis badge first.
    diag_badge <- if (isTRUE(row$is_mactel)) {
      tags$span("MacTel", class = "badge rounded-pill bg-danger me-1")
    } else if (isTRUE(row$is_hsan1)) {
      tags$span("HSAN1", class = "badge rounded-pill bg-warning text-dark me-1")
    } else {
      tags$span("Control", class = "badge rounded-pill bg-secondary me-1")
    }

    # Group-membership badges from the flag columns. Skip HSAN1_variant when
    # the diagnosis badge already shows HSAN1 (HSAN1-only sample).
    tag_cols <- names(SAMPLE_TAG_COLS)
    if (isTRUE(row$is_hsan1) && !isTRUE(row$is_mactel))
      tag_cols <- setdiff(tag_cols, "HSAN1_variant")
    group_badges <- lapply(tag_cols, function(col) {
      if (!(col %in% names(row))) return(NULL)
      val <- suppressWarnings(as.numeric(row[[col]]))
      if (is.na(val) || val != 1) return(NULL)
      tags$span(SAMPLE_TAG_COLS[[col]],
                class = "badge rounded-pill bg-primary me-1")
    })
    group_badges <- Filter(Negate(is.null), group_badges)

    div(class = "mb-2",
        tags$span("Tags: ", class = "text-muted small me-1"),
        diag_badge, group_badges)
  })

  build_sample_dt <- function(d) {
    tbl <- d %>%
      dplyr::transmute(
        Gene = link_gene(SYMBOL), Tier = Tier,
        Variant = link_variant(CHROM, POS, REF, ALT),
        HGVSc, HGVSp = HGVSp_short,
        Impact = IMPACT, Type = TYPE, CADD = round(CADD, 1),
        REVEL = round(REVEL, 3), AlphaMissense = am_class,
        ClinVar = CLNSIG_clean, gnomAD_AF = signif(gnomad_AF, 3),
        Inheritance = inheritance)
    DT::datatable(tbl, rownames = FALSE,
                  selection = "none", escape = FALSE,
                  options = list(pageLength = 15, scrollX = TRUE,
                                 headerCallback = header_tips_cb())) %>%
      DT::formatStyle("ClinVar",
                      backgroundColor = DT::styleEqual(
                        c("Pathogenic", "Pathogenic/Likely_pathogenic",
                          "Likely_pathogenic"),
                        c("#FDDEDE", "#F7D6F7", "#FDE8D8"))) %>%
      DT::formatStyle("AlphaMissense",
                      backgroundColor = DT::styleEqual(
                        "likely_pathogenic", "#FDE8D8")) %>%
      DT::formatStyle("Impact",
                      backgroundColor = DT::styleEqual("HIGH", "#FDDEDE"))
  }

  output$sample_table_priority <- DT::renderDT({
    validate(need(!is.null(input$sample_pick) && input$sample_pick != "",
                  "Choose a sample from the dropdown."))
    d <- sample_data_priority()
    validate(need(nrow(d) > 0,
                  "No variants for this sample pass the current filters."))
    build_sample_dt(d)
  })

  output$sample_table_all <- DT::renderDT({
    validate(need(!is.null(input$sample_pick) && input$sample_pick != "",
                  "Choose a sample from the dropdown."))
    d <- sample_data_all()
    validate(need(nrow(d) > 0, "This sample carries no variants."))
    build_sample_dt(d)
  })

  # ---- click a table row -> detail modal ------------------------------------
  # Row-selection indices from DT refer to the data in its original order, so we
  # re-derive each table's underlying data frame to map row -> variant / gene.
  # Variant tables open a variant-level modal (detail + protein lollipop); the
  # gene-summary table opens the gene-description modal.

  modal_variant <- reactiveVal(NULL)   # one-row data frame of the clicked variant
  modal_gene     <- reactiveVal(NULL)   # gene whose lollipop is on screen

  # ---- shareable bookmark links ---------------------------------------------
  # URL bookmarking (enableBookmarking = "url" on shinyApp) encodes the sidebar
  # filters and active tab into a copy-pasteable link. We additionally save the
  # open gene/variant when the link is generated from inside a protein modal, so
  # the recipient reopens exactly the same view.
  bookmark_with_variant <- reactiveVal(FALSE)

  # Don't pollute the URL with one-shot click/event inputs or the file upload
  # (a local path that wouldn't transfer, and would re-fire on restore).
  setBookmarkExclude(c(
    "upload", "share_link", "share_from_modal", "reset_filters",
    "cell_gene", "cell_variant", "cell_sample",
    "gene_view_variants", "gene_view_lollipop", "open_url"
  ))

  # Sidebar button: link with filters + tab only (no variant).
  observeEvent(input$share_link, {
    bookmark_with_variant(FALSE)
    session$doBookmark()
  })
  # Modal button: link that also reopens the current gene/variant.
  observeEvent(input$share_from_modal, {
    bookmark_with_variant(TRUE)
    session$doBookmark()
  })

  onBookmark(function(state) {
    if (isTRUE(bookmark_with_variant())) {
      state$values$bm_gene <- modal_gene()
      r <- modal_variant()
      state$values$bm_vkey <- if (!is.null(r) && nrow(r) > 0)
        paste(r$CHROM, r$POS, r$REF, r$ALT) else NA_character_
    }
  })

  onBookmarked(function(url) {
    bookmark_with_variant(FALSE)
    showBookmarkUrlModal(url)   # modal with the URL pre-selected for copying
  })

  # On restore the data-load observer above repopulates choices and resets the
  # tier/genes/sample/CADD inputs, clobbering the restored values — so reapply
  # them here, after the session is fully restored, then reopen any saved modal.
  onRestored(function(state) {
    ip <- state$input
    if (!is.null(ip$tier))
      updateCheckboxGroupInput(session, "tier", selected = ip$tier)
    if (!is.null(ip$genes))
      updateSelectizeInput(session, "genes", selected = ip$genes)
    if (!is.null(ip$sample_pick))
      updateSelectizeInput(session, "sample_pick", selected = ip$sample_pick)
    if (!is.null(ip$cadd))
      updateSliderInput(session, "cadd", value = ip$cadd)

    g <- state$values$bm_gene
    if (!is.null(g) && !is.na(g) && nzchar(g)) {
      df <- isolate(raw())
      vk <- state$values$bm_vkey
      if (!is.null(df)) {
        if (!is.null(vk) && !is.na(vk)) {
          hit <- dplyr::filter(df, SYMBOL == g,
                               paste(CHROM, POS, REF, ALT) == vk)
          if (nrow(hit) > 0) { show_variant_modal(hit[1, ]); return() }
        }
        show_gene_lollipop_modal(g)
      }
    }
  })

  # The variant detail block re-renders whenever the selected variant changes,
  # so clicking a different point in the lollipop updates the text at the top.
  output$variant_detail <- renderUI({
    row <- modal_variant(); req(row)
    fld <- function(label, value) {
      if (is.null(value) || length(value) == 0 || is.na(value) || value == "")
        return(NULL)
      tags$span(tags$strong(paste0(label, ": ")), value,
                style = "margin-right:18px;white-space:nowrap;")
    }

    # SpliceAI: max delta score plus the splice event that drives it, with the
    # four component delta scores broken out underneath.
    getn <- function(col) if (col %in% names(row))
      suppressWarnings(as.numeric(row[[col]])) else NA_real_
    ds <- c("acceptor gain" = getn("SpliceAI_pred_DS_AG"),
            "acceptor loss" = getn("SpliceAI_pred_DS_AL"),
            "donor gain"    = getn("SpliceAI_pred_DS_DG"),
            "donor loss"    = getn("SpliceAI_pred_DS_DL"))
    ds_max <- suppressWarnings(max(ds, na.rm = TRUE))
    spliceai_disp <- if (!is.null(row$SpliceAI_max) && !is.na(row$SpliceAI_max)) {
      if (is.finite(ds_max) && ds_max > 0)
        sprintf("%s (%s)", round(row$SpliceAI_max, 3), names(ds)[which.max(ds)])
      else round(row$SpliceAI_max, 3)
    } else NA
    spliceai_breakdown <- if (any(!is.na(ds))) {
      tags$div(class = "text-muted small", style = "margin-top:-4px;",
        sprintf("SpliceAI delta scores — acceptor gain %.2f / acceptor loss %.2f / donor gain %.2f / donor loss %.2f",
                ds[1], ds[2], ds[3], ds[4]))
    } else NULL

    tagList(
      tags$div(
        style = "line-height:1.9;",
        fld("Variant", sprintf("%s:%s %s>%s",
                               row$CHROM, row$POS, row$REF, row$ALT)),
        fld("HGVSc", row$HGVSc),
        fld("HGVSp", row$HGVSp_short),
        tags$br(),
        fld("Impact", as.character(row$IMPACT)),
        fld("Type", as.character(row$TYPE)),
        fld("CADD", if (!is.na(row$CADD)) round(row$CADD, 1) else NA),
        fld("REVEL", if (!is.na(row$REVEL)) round(row$REVEL, 3) else NA),
        fld("AlphaMissense", row$am_class),
        fld("SpliceAI", spliceai_disp),
        tags$br(),
        fld("ClinVar", as.character(row$CLNSIG_clean)),
        fld("gnomAD AF", if (!is.na(row$gnomad_AF)) signif(row$gnomad_AF, 3) else NA),
        fld("Inheritance", row$inheritance),
        spliceai_breakdown
      ),
      tags$hr()
    )
  })

  # External-resource buttons (AlphaFold structure + ClinVar record). Both fire
  # a server-side handler so they open in the real default browser rather than
  # RStudio's blank built-in Viewer. Re-renders with the active variant so the
  # ClinVar link tracks the selected variant.
  output$variant_links <- renderUI({
    row <- modal_variant(); req(row)
    open_btn <- function(url, label, icon) tags$button(
      tagList(bsicons::bs_icon(icon), paste0(" ", label)),
      type = "button", class = "btn btn-sm btn-outline-primary mb-2 me-2",
      onclick = sprintf(
        "Shiny.setInputValue('open_url','%s',{priority:'event'});", .jsesc(url)))

    # AlphaFold: keyed by the gene's UniProt accession from the domain table.
    uni <- if (!is.null(PROTEIN_DOMAINS)) {
      u <- PROTEIN_DOMAINS$UniProt[PROTEIN_DOMAINS$SYMBOL == row$SYMBOL]
      u <- u[!is.na(u) & u != ""]
      if (length(u)) u[1] else NA
    } else NA
    af_btn <- if (!is.na(uni)) open_btn(
      sprintf("https://alphafold.ebi.ac.uk/entry/%s", uni),
      sprintf("AlphaFold structure (%s)", uni), "box") else NULL

    # ClinVar: deep-link to the variation record when a CLNVID is present.
    cv_id <- if ("CLNVID" %in% names(row))
      suppressWarnings(as.character(row$CLNVID)) else NA
    cv_id <- if (length(cv_id) && !is.na(cv_id) && nzchar(cv_id) &&
                 cv_id != "NA") cv_id else NA
    cv_btn <- if (!is.na(cv_id)) open_btn(
      sprintf("https://www.ncbi.nlm.nih.gov/clinvar/variation/%s/", cv_id),
      sprintf("ClinVar record (%s)", cv_id), "clipboard2-pulse") else NULL

    if (is.null(af_btn) && is.null(cv_btn)) return(NULL)
    div(class = "mb-1", af_btn, cv_btn)
  })

  # Samples carrying the active variant, as clickable chips that jump to the
  # sample explorer. Re-renders when the active variant changes.
  output$variant_samples <- renderUI({
    row <- modal_variant(); req(row)
    carriers <- raw() %>%
      dplyr::filter(CHROM == row$CHROM, POS == row$POS,
                    REF == row$REF, ALT == row$ALT) %>%
      dplyr::distinct(family_id) %>%
      dplyr::arrange(family_id) %>%
      dplyr::pull(family_id)
    if (length(carriers) == 0) return(NULL)

    # Diagnosis colour: red MacTel, blue HSAN1, purple both, grey neither.
    diag_colour <- function(fid) {
      if (is.null(SAMPLE_INFO)) return("#9E9E9E")
      r <- SAMPLE_INFO[SAMPLE_INFO$family_id == fid, , drop = FALSE]
      if (nrow(r) == 0) return("#9E9E9E")
      m <- isTRUE(r$is_mactel[1]); h <- isTRUE(r$is_hsan1[1])
      if (m && h) "#9467BD" else if (m) "#D62728" else if (h) "#1F77B4"
      else "#9E9E9E"
    }
    chips <- lapply(carriers, function(fid) {
      tags$a(href = "#", fid,
             class = "badge rounded-pill me-1",
             style = sprintf(
               "background-color:%s;color:#fff;text-decoration:none;cursor:pointer;",
               diag_colour(fid)),
             onclick = sprintf(
               "Shiny.setInputValue('cell_sample','%s',{priority:'event'});return false;",
               .jsesc(fid)))
    })
    legend_dot <- function(col, lab) tags$span(
      tags$span(style = sprintf(
        "display:inline-block;width:10px;height:10px;border-radius:50%%;background-color:%s;margin-right:3px;",
        col)), lab, class = "me-2")
    tagList(
      tags$p(tags$strong(sprintf("Carried by %d sample%s",
                                 length(carriers),
                                 if (length(carriers) == 1) "" else "s")),
             style = "margin-bottom:2px;"),
      div(class = "small text-muted mb-1",
          legend_dot("#D62728", "MacTel"), legend_dot("#1F77B4", "HSAN1"),
          legend_dot("#9467BD", "MacTel + HSAN1"),   legend_dot("#9E9E9E", "Control")),
      div(class = "mb-2", chips)
    )
  })

  show_variant_modal <- function(row) {
    if (is.null(row) || nrow(row) == 0) return(invisible())
    modal_variant(row)
    modal_gene(row$SYMBOL)

    tier <- if (!is.null(row$Tier)) row$Tier else NA

    # gene description (from the bundled gene-info table), shown above the plot
    ginfo <- if (!is.null(GENE_INFO)) GENE_INFO[GENE_INFO$SYMBOL == row$SYMBOL, ] else NULL
    gene_desc <- if (!is.null(ginfo) && nrow(ginfo) > 0 &&
                     !is.na(ginfo$Gene_Description[1]) &&
                     ginfo$Gene_Description[1] != "") {
      tagList(
        tags$p(tags$strong("Gene description"),
               style = "margin-bottom:2px;"),
        tags$p(ginfo$Gene_Description[1], class = "text-muted",
               style = "font-size:0.9rem;")
      )
    } else NULL

    showModal(modalDialog(
      title = tagList(
        tags$span(row$SYMBOL, style = "font-weight:700;font-size:1.2rem;"),
        if (!is.null(tier) && !is.na(tier))
          tags$span(tier, class = "badge bg-secondary",
                    style = "margin-left:8px;vertical-align:middle;")
      ),
      uiOutput("variant_detail"),
      gene_desc,
      uiOutput("variant_links"),
      uiOutput("variant_samples"),
      plotly::plotlyOutput("lollipop", height = 430),
      helpText("Lollipops show every variant in this gene across the loaded ",
               "data (height = CADD, colour = ClinVar, size = number of ",
               "samples). Boxes are Pfam protein domains; the selected ",
               "variant is ringed. Hover a point for detail, or click one to ",
               "switch the summary above to that variant."),
      easyClose = TRUE, size = "xl",
      footer = tagList(
        downloadButton("report_dl",
                       tagList(bsicons::bs_icon("file-earmark-arrow-down"),
                               " Download report"),
                       class = "btn btn-primary", icon = NULL),
        tags$button(
          tagList(bsicons::bs_icon("link-45deg"), " Copy share link"),
          class = "btn btn-outline-primary",
          title = "Get a link that reopens this gene/variant with the current filters.",
          onclick = paste0("Shiny.setInputValue('share_from_modal', Math.random(),",
                           " {priority:'event'}); return false;")),
        modalButton("Close"))
    ))
  }

  # Gene-level landing for the protein lollipop: shows the plot for every variant
  # in the gene with no variant pre-selected (the detail header stays blank until
  # a point is clicked). Same modal/visuals as the variant view, minus the ring.
  show_gene_lollipop_modal <- function(symbol) {
    if (is.null(symbol) || is.na(symbol) || symbol == "") return(invisible())
    modal_variant(NULL)
    modal_gene(symbol)

    info <- if (!is.null(GENE_INFO)) GENE_INFO[GENE_INFO$SYMBOL == symbol, ] else NULL
    tier <- if (!is.null(info) && nrow(info) > 0) info$Tier[1] else NA
    gene_desc <- if (!is.null(info) && nrow(info) > 0 &&
                     !is.na(info$Gene_Description[1]) &&
                     info$Gene_Description[1] != "") {
      tagList(
        tags$p(tags$strong("Gene description"), style = "margin-bottom:2px;"),
        tags$p(info$Gene_Description[1], class = "text-muted",
               style = "font-size:0.9rem;")
      )
    } else NULL

    showModal(modalDialog(
      title = tagList(
        tags$span(symbol, style = "font-weight:700;font-size:1.2rem;"),
        if (!is.null(tier) && !is.na(tier))
          tags$span(tier, class = "badge bg-secondary",
                    style = "margin-left:8px;vertical-align:middle;")
      ),
      uiOutput("variant_detail"),
      gene_desc,
      uiOutput("variant_links"),
      uiOutput("variant_samples"),
      plotly::plotlyOutput("lollipop", height = 430),
      helpText("Lollipops show every variant in this gene across the loaded ",
               "data (height = CADD, colour = ClinVar, size = number of ",
               "samples). Boxes are Pfam protein domains. Hover a point for ",
               "detail, or click one to load that variant above."),
      easyClose = TRUE, size = "xl",
      footer = tagList(
        downloadButton("report_dl",
                       tagList(bsicons::bs_icon("file-earmark-arrow-down"),
                               " Download report"),
                       class = "btn btn-primary", icon = NULL),
        tags$button(
          tagList(bsicons::bs_icon("link-45deg"), " Copy share link"),
          class = "btn btn-outline-primary",
          title = "Get a link that reopens this gene/variant with the current filters.",
          onclick = paste0("Shiny.setInputValue('share_from_modal', Math.random(),",
                           " {priority:'event'}); return false;")),
        modalButton("Close"))
    ))
  }

  # ---- one-click report (self-contained HTML) ------------------------------
  # Build a standalone, emailable HTML summary for the gene/variant currently
  # shown in the protein modal. The lollipop is embedded as a base64 PNG so the
  # file is fully self-contained — no internet or extra files needed to view it.
  make_report_html <- function(gene, row) {
    esc <- function(x) htmltools::htmlEscape(as.character(x))
    has_row <- !is.null(row) && nrow(row) > 0

    ginfo <- if (!is.null(GENE_INFO)) GENE_INFO[GENE_INFO$SYMBOL == gene, ] else NULL
    tier  <- if (!is.null(ginfo) && nrow(ginfo) > 0) ginfo$Tier[1]
             else if (has_row && !is.null(row$Tier)) row$Tier else NA
    gdesc <- if (!is.null(ginfo) && nrow(ginfo) > 0 &&
                 !is.na(ginfo$Gene_Description[1]) &&
                 nzchar(ginfo$Gene_Description[1])) ginfo$Gene_Description[1] else NULL

    # Lollipop -> base64-embedded PNG.
    gdf <- dplyr::filter(raw(), SYMBOL == gene)
    ddf <- if (!is.null(PROTEIN_DOMAINS))
      dplyr::filter(PROTEIN_DOMAINS, SYMBOL == gene) else NULL
    sel_key <- if (has_row) paste(row$CHROM, row$POS, row$REF, row$ALT) else NULL
    p <- tryCatch(plot_variant_lollipop(gdf, ddf, gene, sel_key),
                  error = function(e) NULL)
    plot_html <- "<p class='muted'>No protein-coding positions to plot for this gene.</p>"
    if (!is.null(p)) {
      tmp <- tempfile(fileext = ".png")
      ok <- tryCatch({
        ggplot2::ggsave(tmp, p, width = 9, height = 3.8, dpi = 130, bg = "white")
        TRUE
      }, error = function(e) FALSE)
      if (ok && file.exists(tmp)) {
        uri <- paste0("data:image/png;base64,",
                      jsonlite::base64_enc(readBin(tmp, "raw", file.info(tmp)$size)))
        plot_html <- sprintf('<img class="plot" src="%s" alt="Protein lollipop"/>', uri)
        unlink(tmp)
      }
    }

    if (has_row) {
      getn <- function(col) if (col %in% names(row))
        suppressWarnings(as.numeric(row[[col]])) else NA_real_
      ds <- c("acceptor gain" = getn("SpliceAI_pred_DS_AG"),
              "acceptor loss" = getn("SpliceAI_pred_DS_AL"),
              "donor gain"    = getn("SpliceAI_pred_DS_DG"),
              "donor loss"    = getn("SpliceAI_pred_DS_DL"))
      spliceai_disp <- if (!is.null(row$SpliceAI_max) && !is.na(row$SpliceAI_max)) {
        m <- suppressWarnings(max(ds, na.rm = TRUE))
        if (is.finite(m) && m > 0)
          sprintf("%s (%s)", round(row$SpliceAI_max, 3), names(ds)[which.max(ds)])
        else as.character(round(row$SpliceAI_max, 3))
      } else NA
      fields <- list(
        c("Variant",       sprintf("%s:%s %s>%s", row$CHROM, row$POS, row$REF, row$ALT)),
        c("HGVSc",         row$HGVSc),
        c("HGVSp",         row$HGVSp_short),
        c("Impact",        as.character(row$IMPACT)),
        c("Type",          as.character(row$TYPE)),
        c("CADD",          if (!is.na(row$CADD)) round(row$CADD, 1) else NA),
        c("REVEL",         if (!is.na(row$REVEL)) round(row$REVEL, 3) else NA),
        c("AlphaMissense", row$am_class),
        c("SpliceAI",      spliceai_disp),
        c("ClinVar",       as.character(row$CLNSIG_clean)),
        c("gnomAD AF",     if (!is.na(row$gnomad_AF)) signif(row$gnomad_AF, 3) else NA),
        c("Inheritance",   row$inheritance)
      )
      kv <- vapply(fields, function(f) {
        v <- f[[2]]
        if (is.null(v) || length(v) == 0 || is.na(v) || v == "") return("")
        sprintf("<tr><th>%s</th><td>%s</td></tr>", esc(f[[1]]), esc(v))
      }, character(1))
      body_html <- sprintf("<table class='kv'>%s</table>", paste(kv, collapse = ""))

      carriers <- gdf %>%
        dplyr::filter(CHROM == row$CHROM, POS == row$POS,
                      REF == row$REF, ALT == row$ALT) %>%
        dplyr::distinct(family_id) %>% dplyr::arrange(family_id) %>%
        dplyr::pull(family_id)
      diag_colour <- function(fid) {
        if (is.null(SAMPLE_INFO)) return("#9E9E9E")
        r <- SAMPLE_INFO[SAMPLE_INFO$family_id == fid, , drop = FALSE]
        if (nrow(r) == 0) return("#9E9E9E")
        m <- isTRUE(r$is_mactel[1]); h <- isTRUE(r$is_hsan1[1])
        if (m && h) "#9467BD" else if (m) "#D62728" else if (h) "#1F77B4" else "#9E9E9E"
      }
      if (length(carriers) > 0) {
        pills <- paste(vapply(carriers, function(fid)
          sprintf('<span class="pill" style="background:%s">%s</span>',
                  diag_colour(fid), esc(fid)), character(1)), collapse = "")
        body_html <- paste0(body_html,
          sprintf("<h3>Carried by %d sample%s</h3>", length(carriers),
                  if (length(carriers) == 1) "" else "s"),
          "<p class='legend'>",
          "<span class='dot' style='background:#D62728'></span>MacTel",
          "<span class='dot' style='background:#1F77B4'></span>HSAN1",
          "<span class='dot' style='background:#9467BD'></span>MacTel + HSAN1",
          "<span class='dot' style='background:#9E9E9E'></span>Control</p>",
          sprintf("<div class='pills'>%s</div>", pills))
      }
      title_line <- sprintf("%s &mdash; %s:%s %s&gt;%s", esc(gene), esc(row$CHROM),
                            esc(row$POS), esc(row$REF), esc(row$ALT))
    } else {
      vsum <- gdf %>%
        dplyr::group_by(CHROM, POS, REF, ALT) %>%
        dplyr::summarise(HGVSp   = dplyr::first(HGVSp_short),
                         Impact  = as.character(dplyr::first(IMPACT)),
                         CADD    = suppressWarnings(max(CADD, na.rm = TRUE)),
                         ClinVar = as.character(dplyr::first(CLNSIG_clean)),
                         Samples = dplyr::n_distinct(family_id),
                         .groups = "drop") %>%
        dplyr::arrange(dplyr::desc(CADD))
      vr <- vapply(seq_len(nrow(vsum)), function(i) {
        r <- vsum[i, ]
        sprintf("<tr><td>%s:%s %s&gt;%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td><td>%s</td></tr>",
                esc(r$CHROM), esc(r$POS), esc(r$REF), esc(r$ALT),
                esc(r$HGVSp), esc(r$Impact),
                esc(if (is.finite(r$CADD)) round(r$CADD, 1) else ""),
                esc(r$ClinVar), esc(r$Samples))
      }, character(1))
      body_html <- sprintf(
        "<h3>%d variant%s in this gene</h3><table class='vt'><thead><tr><th>Variant</th><th>HGVSp</th><th>Impact</th><th>CADD</th><th>ClinVar</th><th>Samples</th></tr></thead><tbody>%s</tbody></table>",
        nrow(vsum), if (nrow(vsum) == 1) "" else "s", paste(vr, collapse = ""))
      title_line <- esc(gene)
    }

    tier_badge <- if (!is.null(tier) && !is.na(tier))
      sprintf("<span class='tier'>%s</span>", esc(tier)) else ""
    desc_html <- if (!is.null(gdesc))
      sprintf("<h3>Gene description</h3><p>%s</p>", esc(gdesc)) else ""

    paste0(
      "<!DOCTYPE html><html lang='en'><head><meta charset='utf-8'>",
      sprintf("<title>MacTel report — %s</title>", esc(gene)),
      "<style>",
      "body{font-family:-apple-system,Segoe UI,Roboto,Helvetica,Arial,sans-serif;color:#222;max-width:900px;margin:24px auto;padding:0 18px;line-height:1.5;}",
      "h1{font-size:1.5rem;margin:0;color:#1F4E79;}",
      "h3{font-size:1.05rem;margin:18px 0 6px;color:#1F4E79;border-bottom:1px solid #eee;padding-bottom:3px;}",
      ".sub{color:#666;font-size:0.85rem;margin:2px 0 14px;}",
      ".tier{background:#1F4E79;color:#fff;border-radius:10px;padding:2px 8px;font-size:0.8rem;margin-left:8px;vertical-align:middle;}",
      "table.kv{border-collapse:collapse;}table.kv th{text-align:left;padding:3px 14px 3px 0;color:#555;font-weight:600;white-space:nowrap;vertical-align:top;}table.kv td{padding:3px 0;}",
      "table.vt{border-collapse:collapse;width:100%;font-size:0.88rem;}table.vt th,table.vt td{border:1px solid #e3e3e3;padding:4px 8px;text-align:left;}table.vt th{background:#f4f7fa;}",
      "img.plot{max-width:100%;height:auto;border:1px solid #eee;border-radius:6px;margin-top:6px;}",
      ".muted{color:#888;}",
      ".pills{margin-top:4px;}.pill{display:inline-block;color:#fff;border-radius:10px;padding:2px 9px;margin:0 4px 4px 0;font-size:0.8rem;}",
      ".legend{font-size:0.8rem;color:#666;margin:4px 0;}.dot{display:inline-block;width:10px;height:10px;border-radius:50%;margin:0 3px 0 10px;}",
      ".foot{margin-top:24px;color:#999;font-size:0.78rem;border-top:1px solid #eee;padding-top:8px;}",
      "</style></head><body>",
      sprintf("<h1>%s%s</h1>", title_line, tier_badge),
      sprintf("<p class='sub'>MacTel Variant Explorer report &middot; generated %s</p>", Sys.Date()),
      body_html,
      desc_html,
      "<h3>Protein lollipop</h3>", plot_html,
      "<p class='foot'>Generated by the MacTel Variant Explorer. Scores are ",
      "computational predictions &mdash; interpret in clinical context.</p>",
      "</body></html>"
    )
  }

  output$report_dl <- downloadHandler(
    filename = function() {
      g <- modal_gene(); r <- modal_variant()
      g <- if (is.null(g)) "gene" else g
      if (!is.null(r) && nrow(r) > 0)
        sprintf("MacTel_report_%s_%s-%s_%s.html", g, r$CHROM, r$POS, Sys.Date())
      else
        sprintf("MacTel_report_%s_%s.html", g, Sys.Date())
    },
    content = function(file) {
      g <- modal_gene()
      if (is.null(g)) {
        writeLines("<html><body><p>No gene selected.</p></body></html>", file)
        return(invisible())
      }
      writeLines(make_report_html(g, modal_variant()), file)
    }
  )

  output$lollipop <- plotly::renderPlotly({
    gene <- modal_gene(); req(gene)
    gdf <- dplyr::filter(raw(), SYMBOL == gene)
    ddf <- if (!is.null(PROTEIN_DOMAINS))
      dplyr::filter(PROTEIN_DOMAINS, SYMBOL == gene) else NULL
    # Only ring a variant once one has been selected (blank on the gene landing).
    row <- modal_variant()
    sel_key <- if (!is.null(row) && nrow(row) > 0)
      paste(row$CHROM, row$POS, row$REF, row$ALT) else NULL
    p <- plot_variant_lollipop(gdf, ddf, gene, sel_key)
    validate(need(!is.null(p),
                  "No protein-coding (amino-acid) positions to plot for this gene."))
    plotly::ggplotly(p, tooltip = "text", source = "lollipop") %>%
      plotly::event_register("plotly_click")
  })

  # Clicking a point in the lollipop switches the active variant.
  observeEvent(plotly::event_data("plotly_click", source = "lollipop"), {
    ed   <- plotly::event_data("plotly_click", source = "lollipop")
    gene <- modal_gene()
    req(ed, gene)
    key <- ed$key
    if (is.null(key) || length(key) == 0 || is.na(key[[1]])) return()
    hit <- raw() %>%
      dplyr::filter(SYMBOL == gene,
                    paste(CHROM, POS, REF, ALT) == key[[1]])
    if (nrow(hit) > 0) modal_variant(hit[1, ])
  })

  # Clicking a Gene cell -> gene-description modal.
  observeEvent(input$cell_gene, {
    show_gene_modal(input$cell_gene)
  })

  # "View all variants" in the gene modal -> set the gene filter to that gene
  # (exactly as if picked from the sidebar) and jump to the Variant table.
  # Global filters still apply, since this only adds to the gene selection.
  observeEvent(input$gene_view_variants, {
    removeModal()
    df <- raw(); req(df)
    updateSelectizeInput(session, "genes",
                         choices = sort(unique(df$SYMBOL)),
                         selected = input$gene_view_variants,
                         server = TRUE)
    bslib::nav_select("main_tabs", "Variant table")
  })

  # "View variants on protein" in the gene modal -> open the protein lollipop for
  # the whole gene with no variant pre-selected (detail header blank until click).
  observeEvent(input$gene_view_lollipop, {
    show_gene_lollipop_modal(input$gene_view_lollipop)
  })

  # Clicking a Variant cell -> variant detail + protein lollipop modal.
  observeEvent(input$cell_variant, {
    parts <- strsplit(input$cell_variant, "\\|\\|")[[1]]
    if (length(parts) != 4) return()
    hit <- raw() %>%
      dplyr::filter(CHROM == parts[1], as.character(POS) == parts[2],
                    REF == parts[3], ALT == parts[4])
    if (nrow(hit) > 0) show_variant_modal(hit[1, ])
  })

  # Clicking a Sample cell (table or variant-modal chip) -> close any modal and
  # jump to the Sample explorer for that sample.
  observeEvent(input$cell_sample, {
    removeModal()
    updateSelectizeInput(session, "sample_pick",
                         selected = input$cell_sample)
    bslib::nav_select("main_tabs", "Sample explorer")
  })

  # Open a URL (AlphaFold, ClinVar, …) in the OS default browser, bypassing
  # RStudio's browser option (which otherwise traps it in its blank Viewer).
  observeEvent(input$open_url, {
    url <- input$open_url
    if (is.null(url) || !nzchar(url)) return()
    os <- Sys.info()[["sysname"]]
    tryCatch(
      if (os == "Darwin") {
        system2("open", shQuote(url), wait = FALSE)
      } else if (os == "Windows") {
        system2("cmd", c("/c", "start", '""', shQuote(url)), wait = FALSE)
      } else {
        system2("xdg-open", shQuote(url), wait = FALSE)
      },
      error = function(e) utils::browseURL(url)
    )
  })

  # ---- downloads ------------------------------------------------------------
  output$dl_table <- downloadHandler(
    filename = function() sprintf("filtered_variants_%s.csv", Sys.Date()),
    content  = function(file) readr::write_csv(display_cols(filtered(), links = FALSE), file)
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
      d <- sample_data_all() %>%
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

shinyApp(ui, server, enableBookmarking = "url")

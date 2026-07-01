# =============================================================================
# MacTel Variant Explorer — Shiny app
#
# Interactive exploration of Cavalier candidate-variant output.
# Expects the 71-column Cavalier CSV structure (see R/load_data.R).
#
# Run locally:
#   shiny::runApp("path/to/mactel_variant_explorer")
# or open app.R in RStudio and click "Run App".
#
# Data source priority:
#   1. data/candidate_variants.csv  (real data, gitignored — used if present)
#   2. a file uploaded via the sidebar
# No data ships with the repo; point the app at your own Cavalier CSV.
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
library(shinyFiles)
library(patchwork)   # attach: the | and / layout operators need it on the path

# Files in R/ are auto-sourced by Shiny, but source explicitly so the app also
# works when launched via shiny::runApp() from any working directory.
app_dir <- tryCatch(dirname(sys.frame(1)$ofile), error = function(e) ".")
for (f in list.files(file.path(app_dir, "R"), pattern = "\\.R$", full.names = TRUE)) {
  source(f, local = FALSE)
}

DEFAULT_REAL    <- file.path(app_dir, "data", "candidate_variants_tier1_2.csv")
DEFAULT_EXAMPLE <- file.path(app_dir, "data", "example_variants.csv")
startup_path    <- if (file.exists(DEFAULT_REAL)) DEFAULT_REAL else DEFAULT_EXAMPLE

# DEBUG mode -----------------------------------------------------------------
# When OFF (the default) the app starts with NO variant file loaded — the user
# must explicitly upload/select their own Cavalier CSV. This guards against
# accidentally working with a stale or bundled file.
# When ON, the app auto-loads `startup_path` at launch (handy for development).
# Turn it on with either:
#   options(mactel.debug = TRUE)         # in R, before runApp()
#   MACTEL_DEBUG=TRUE  (environment variable)
DEBUG <- isTRUE(getOption("mactel.debug", FALSE)) ||
  identical(toupper(Sys.getenv("MACTEL_DEBUG", "")), "TRUE")

# Optional folder of per-sample IGV reports (one sub-folder per individual,
# each holding <SAMPLE>.igv_report.html). Auto-detected if bundled under
# data/igv_reports; otherwise the user points the app at it from the sidebar.
# Holds real sample data, so it is never committed — pointed-to at runtime.
DEFAULT_IGV_DIR <- file.path(app_dir, "data", "igv_reports")
igv_startup     <- if (dir.exists(DEFAULT_IGV_DIR)) DEFAULT_IGV_DIR else ""

# Gene info (Tier + descriptions; gene symbols only, no participant data) ------
# The curated candidate-gene list is versioned: each release is a date-stamped
# TSV (gene_info_YYYY-MM-DD.tsv) under data/gene_lists/. The user picks a version
# from the sidebar; the newest is the default. Falls back to a single legacy
# data/gene_info.tsv if the folder is absent.
GENE_LIST_DIR <- file.path(app_dir, "data", "gene_lists")

# Returns a named character vector: display label (version date, newest first)
# -> file path. Labels are the YYYY-MM-DD pulled from the filename when present.
list_gene_lists <- function() {
  files <- if (dir.exists(GENE_LIST_DIR))
    sort(list.files(GENE_LIST_DIR, pattern = "\\.tsv$", full.names = TRUE),
         decreasing = TRUE) else character(0)
  if (length(files) == 0) {
    legacy <- file.path(app_dir, "data", "gene_info.tsv")
    if (file.exists(legacy)) files <- legacy
  }
  if (length(files) == 0) return(character(0))
  labels <- vapply(files, function(f) {
    d <- regmatches(basename(f), regexpr("[0-9]{4}-[0-9]{2}-[0-9]{2}", basename(f)))
    if (length(d) == 1 && nzchar(d)) d else tools::file_path_sans_ext(basename(f))
  }, character(1))
  stats::setNames(files, labels)
}

GENE_LISTS     <- list_gene_lists()
GENE_INFO_PATH <- if (length(GENE_LISTS)) unname(GENE_LISTS[1]) else
                  file.path(app_dir, "data", "gene_info.tsv")
GENE_INFO      <- load_gene_info(GENE_INFO_PATH)

# Pfam protein domains, pre-fetched offline (gene symbols only, no patient data).
PROTEIN_DOMAINS <- load_protein_domains(file.path(app_dir, "data",
                                                  "protein_domains.tsv"))

# The embedded 3D structure viewer needs the optional r3dmol package. When it is
# absent the app still runs — the structure section is simply omitted.
HAS_R3DMOL <- requireNamespace("r3dmol", quietly = TRUE)

# Per-sample info (case/control status + data-group flags). Prefer the real
# sheet when present (git-ignored); otherwise the de-identified example.
SAMPLE_INFO_REAL    <- file.path(app_dir, "data", "all_samples_fixed.txt")
SAMPLE_INFO_EXAMPLE <- file.path(app_dir, "data", "example_sample_info.tsv")
SAMPLE_INFO <- load_sample_info(
  if (file.exists(SAMPLE_INFO_REAL)) SAMPLE_INFO_REAL else SAMPLE_INFO_EXAMPLE)

# family_id -> alternate-ID lookups for the sample-ID display toggle. Built once
# from SAMPLE_INFO. The variant data and every internal join stay keyed on
# family_id; only the *displayed* label changes. Falls back to family_id when a
# mapping is missing or blank.
AID_LOOKUP <- character(0)
PID_LOOKUP <- character(0)
if (!is.null(SAMPLE_INFO) && "family_id" %in% names(SAMPLE_INFO)) {
  if ("AID" %in% names(SAMPLE_INFO)) {
    AID_LOOKUP <- stats::setNames(as.character(SAMPLE_INFO$AID),
                                  SAMPLE_INFO$family_id)
  }
  if ("Patient_ID" %in% names(SAMPLE_INFO)) {
    PID_LOOKUP <- stats::setNames(as.character(SAMPLE_INFO$Patient_ID),
                                  SAMPLE_INFO$family_id)
  }
}

# family_id -> diagnosis group, for the stacked top-genes bar chart. One group
# per sample, mirroring diag_colour(): both -> "MacTel + HSAN1", else MacTel /
# HSAN1 / Control. NULL when no sample info is available.
DIAG_GROUP_LOOKUP <- NULL
if (!is.null(SAMPLE_INFO) &&
    all(c("family_id", "is_mactel", "is_hsan1") %in% names(SAMPLE_INFO))) {
  grp <- ifelse(SAMPLE_INFO$is_mactel & SAMPLE_INFO$is_hsan1, "MacTel + HSAN1",
         ifelse(SAMPLE_INFO$is_mactel, "MacTel",
         ifelse(SAMPLE_INFO$is_hsan1,  "HSAN1", "Control")))
  DIAG_GROUP_LOOKUP <- stats::setNames(grp, SAMPLE_INFO$family_id)
}

# Map a vector of family_ids to display labels for the chosen format ("AID" or
# "Patient ID"). Vectorised; preserves order; blanks/NAs fall back to family_id.
format_sample_id <- function(fid, fmt = "AID") {
  fid <- as.character(fid)
  lk  <- if (identical(fmt, "Patient ID")) PID_LOOKUP else AID_LOOKUP
  if (length(lk) == 0) return(fid)
  out <- unname(lk[fid])
  ifelse(is.na(out) | !nzchar(out), fid, out)
}

# --- Excel export ------------------------------------------------------------
# Sensible default column selection for the "Export to Excel" picker. The picker
# offers every column present in the loaded data (all original Cavalier fields
# plus the app's derived ones); these are just the ones ticked by default. Only
# columns that actually exist are used, so this is safe across CSV variants.
# These mirror the columns shown in the gene report / gene-landing table
# (Sample, Variant, HGVSp protein change without the transcript prefix, Impact,
# CADD, ClinVar class) plus the gene SYMBOL and the variant Consequence.
EXPORT_DEFAULT_COLS <- c(
  "SYMBOL", "family_id", "Variant", "HGVSp_short", "Consequence",
  "IMPACT", "CADD", "CLNSIG_clean")
# Extra defaults when exporting the Priority variants tab.
EXPORT_PRIORITY_EXTRA <- c("n_flags", "why_prioritised",
                           "flag_clinvar", "flag_high", "flag_cadd")

# Write a data frame to a nicely formatted .xlsx: a banded Excel table with a
# styled header row, auto-filter, frozen header + first column, and column
# widths sized to content (capped so long text columns stay readable).
write_variants_xlsx <- function(df, file, sheet = "Variants") {
  wb <- openxlsx::createWorkbook()
  openxlsx::addWorksheet(wb, sheet)
  openxlsx::writeDataTable(wb, sheet, df, tableStyle = "TableStyleMedium2",
                           withFilter = TRUE, bandedRows = TRUE)
  openxlsx::freezePane(wb, sheet, firstActiveRow = 2, firstActiveCol = 2)
  widths <- vapply(seq_along(df), function(i) {
    vals <- as.character(df[[i]]); vals <- vals[!is.na(vals)]
    body <- if (length(vals)) max(nchar(vals)) else 0
    min(max(body, nchar(names(df)[i]), 8) + 2, 45)
  }, numeric(1))
  openxlsx::setColWidths(wb, sheet, cols = seq_along(df), widths = widths)
  openxlsx::saveWorkbook(wb, file, overwrite = TRUE)
}

# Family grouping (from the manifest's Family_ID column, derived from the
# clinical pedigree). FAMILY_OF maps an individual family_id -> its Family_ID;
# FAMILY_MEMBERS maps a Family_ID -> the vector of member family_ids. Only
# samples in a multi-member family carry a Family_ID; singletons are absent.
FAMILY_OF      <- character(0)
FAMILY_MEMBERS <- list()
if (!is.null(SAMPLE_INFO) &&
    all(c("family_id", "Family_ID") %in% names(SAMPLE_INFO))) {
  fam <- SAMPLE_INFO[!is.na(SAMPLE_INFO$Family_ID) & nzchar(SAMPLE_INFO$Family_ID),
                     c("family_id", "Family_ID"), drop = FALSE]
  if (nrow(fam)) {
    FAMILY_OF      <- stats::setNames(as.character(fam$Family_ID),
                                      as.character(fam$family_id))
    FAMILY_MEMBERS <- split(as.character(fam$family_id),
                            as.character(fam$Family_ID))
  }
}

# Family labels for the Family explorer dropdown, ordered by their numeric
# suffix (FAMILY1, FAMILY2, … FAMILY10) rather than lexicographically.
FAMILY_CHOICES <- names(FAMILY_MEMBERS)
if (length(FAMILY_CHOICES)) {
  FAMILY_CHOICES <- FAMILY_CHOICES[
    order(suppressWarnings(as.integer(gsub("\\D", "", FAMILY_CHOICES))),
          FAMILY_CHOICES)]
}

# Family_ID for a single individual family_id, or NA when it has no family.
family_of <- function(fid) {
  fid <- as.character(fid)
  if (length(FAMILY_OF) == 0 || length(fid) != 1 || is.na(fid)) return(NA_character_)
  out <- unname(FAMILY_OF[fid])
  if (is.na(out) || !nzchar(out)) NA_character_ else out
}

# Tier lookup: prefer the richer gene_info table; fall back to gene_tiers.tsv.
TIER_PATH <- file.path(app_dir, "data", "gene_tiers.tsv")
TIER_DF   <- if (!is.null(GENE_INFO)) {
  dplyr::distinct(dplyr::select(GENE_INFO, SYMBOL, Tier))
} else {
  load_gene_tiers(TIER_PATH)
}

# Count curated genes in a tier (for the landing-page summary box). Works on
# whichever gene-list version is active (a tier lookup tibble of SYMBOL + Tier).
tier_gene_count <- function(tier_df, label) {
  if (is.null(tier_df) || !"Tier" %in% names(tier_df)) return(0L)
  sum(tier_df$Tier == label, na.rm = TRUE)
}

# Build a modalDialog describing a single gene from the active gene-info table.
show_gene_modal <- function(symbol, gene_info = GENE_INFO) {
  if (is.null(symbol) || is.na(symbol) || symbol == "") return(invisible())
  info <- if (!is.null(gene_info)) gene_info[gene_info$SYMBOL == symbol, ] else NULL

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

# Load + clean + annotate tier in one step. Pass the tier lookup for the gene
# list version currently selected (defaults to the startup TIER_DF).
load_annotated <- function(path, tier_df = TIER_DF) {
  annotate_tier(load_variants(path), tier_df)
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
  "Only in"       = "Which of the two compared files this variant is unique to",
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

    /* Tab bar — larger, more readable tab labels. */
    .nav-tabs .nav-link {
      padding: 0.4rem 0.8rem !important;
      font-size: 1rem !important;
    }

    /* Header summary stat cards: keep the live count text white. */
    .bslib-grid > div .shiny-text-output { color: inherit !important; }

    /* Larger card sub-headings across the app (e.g. How to use this app,
       The tabs, MacTel gene tiers, and each chart's title). */
    .card-header { font-size: 1.15rem !important; }
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

# One tier "card" for the landing-page tier summary box: a big gene count, the
# tier name, and a short plain-English description of the evidence strength.
tier_box <- function(label, n, desc, accent) div(
  class = "d-flex align-items-start p-3 rounded-3 h-100",
  style = paste0("background:", accent$bg, ";border-left:5px solid ", accent$bar, ";"),
  div(class = "flex-shrink-0 text-center me-3",
      div(style = paste0("font-size:2rem;font-weight:700;line-height:1;color:",
                         accent$bar, ";"),
          n),
      div(class = "small text-body-secondary", "genes")),
  div(
    tags$div(class = "fw-semibold", label),
    tags$div(class = "small text-body-secondary", desc))
)

# A header summary stat: a coloured tile with the icon glued to the left of the
# label + live count. Built with flexbox (not bslib value_box) so the icon
# always sits right next to the text, at any window width.
stat_card <- function(value_id, label, icon, bg) div(
  class = "d-flex align-items-center p-3 rounded-3 h-100 text-white",
  style = paste0("background:", bg, ";"),
  div(class = "flex-shrink-0 d-flex align-items-center justify-content-center me-3",
      style = paste0("width:3rem;height:3rem;border-radius:0.6rem;",
                     "background:rgba(255,255,255,0.22);font-size:1.6rem;"),
      if (inherits(icon, c("html", "shiny.tag", "shiny.tag.list")))
        icon else bsicons::bs_icon(icon)),
  div(class = "lh-1",
      tags$div(label, style = "font-size:0.9rem;font-weight:600;opacity:0.92;"),
      tags$div(textOutput(value_id, inline = TRUE),
               style = "font-size:1.8rem;font-weight:700;margin-top:0.15rem;"))
)

ui <- function(request) page_sidebar(
  title = tags$span(
    class = "d-inline-flex align-items-center",
    dna_icon(40),
    tags$span("MacTel Variant Explorer", class = "ms-2",
              style = "font-size:2.2rem;font-weight:400;letter-spacing:0.2px;")
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
                  accept = c(".csv"), buttonLabel = "Browse…"),
        tags$label("IGV reports folder", class = "control-label d-block"),
        shinyFiles::shinyDirButton(
          "igv_dir_btn",
          label = "Choose folder…",
          title = "Select the folder of per-sample IGV reports",
          icon = bsicons::bs_icon("folder2-open"),
          class = "btn-outline-primary btn-sm w-100"),
        div(class = "small text-muted mt-1 text-break",
            textOutput("igv_dir_label", inline = TRUE)),
        tags$hr(class = "my-2"),
        selectInput("gene_list_pick", "Gene list version",
                    choices = names(GENE_LISTS),
                    selected = if (length(GENE_LISTS)) names(GENE_LISTS)[1] else NULL,
                    width = "100%"),
        tags$hr(class = "my-2"),
        radioButtons("id_format", "Sample ID format",
                     choices = c("AID", "Patient ID"),
                     selected = "AID", inline = TRUE)
      ),

      accordion_panel(
        "Core filters", icon = bsicons::bs_icon("funnel"),
        checkboxGroupInput("sample_group", "Sample group",
                           choices = c("MacTel", "HSAN1", "Controls"),
                           selected = "MacTel", inline = TRUE),
        selectizeInput("exclude_samples", "Exclude samples",
                       choices = NULL, multiple = TRUE,
                       options = list(placeholder = "None excluded")),
        helpText(class = "small text-muted mt-n2",
                 "Drop specific samples from every view (e.g. suspected ",
                 "bad-data samples), regardless of the group selected above."),
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
        tags$label("AlphaMissense", class = "control-label"),
        checkboxInput("exclude_am_benign",
                      "Exclude benign / likely-benign", FALSE),
        sliderInput("cadd", "CADD ≥",
                    min = 0, max = 60, value = 0, step = 1),
        sliderInput("revel", "REVEL ≥",
                    min = 0, max = 1, value = 0, step = 0.05),
        numericInput("gnomad", "gnomAD AF ≤",
                     value = NA, min = 0, max = 1, step = 0.001),
        helpText(class = "small text-muted mt-n2",
                 "Max population allele frequency (0–1). Leave blank for no limit.")
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

    tags$hr(class = "my-2"),
    fileInput("compare_upload", "Compare variants",
              accept = c(".csv"), buttonLabel = "Browse…"),
    helpText(class = "small",
             "Upload a second Cavalier CSV to open a 'Variant comparisons' ",
             "tab listing variants present in only one of the two files."),

    actionButton("filter_share",
                 tagList(bsicons::bs_icon("sliders"), " Share / save filters"),
                 class = "btn-outline-primary btn-sm w-100 mt-3",
                 title = paste("Copy a short code of the current filter settings to",
                               "share with a colleague or save for later — paste a",
                               "code to apply those filters to your own data."))
  ),

  # Header summary stats (custom flexbox cards — see stat_card()).
  layout_columns(
    fill = FALSE,
    col_widths = c(3, 3, 3, 3),
    stat_card("vb_variants", "Variants", "file-earmark-text", "#1F4E79"),
    stat_card("vb_genes",    "Genes",    dna_icon(26),        "#5a7184"),
    stat_card("vb_samples",  "Samples",  "people",            "#3498DB"),
    stat_card("vb_plp",      "ClinVar P/LP", "exclamation-triangle", "#C0392B")
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
              "MacTel study")
          )
        )
      ),
      layout_columns(
        col_widths = c(7, 5),
        card(
          card_header(bsicons::bs_icon("info-circle"),
                      tags$strong(" How to use this app")),
          tags$p(class = "text-body-secondary small mb-3",
            HTML("Each row in the data is a <strong>variant</strong> seen in a <strong>sample</strong>, in a MacTel gene.")),
          landing_step(1, "Filter on the left",
            "Narrow by gene, sample group, predicted severity, or how rare the variant is. The counters up top update live."),
          landing_step(2, "Browse the tabs",
            "Each tab shows the filtered variants a different way — charts, tables, and plots."),
          landing_step(3, "Click anything blue",
            "Genes, variants, and samples are links that open detail views."),
          landing_step(4, "Save or share your filters",
            tagList("The ", tags$strong("Share / save filters"), " button (bottom of ",
                    "the sidebar) turns your current filters into a short code. ",
                    "Send it to a colleague or paste it back later to restore the ",
                    "same view — it works on any loaded dataset.")),
          landing_step(5, "Export the results",
            tagList("Download any table as ", tags$strong("CSV"), ", or open a gene ",
                    "or variant and generate a one-click, self-contained ",
                    tags$strong("HTML report"), " to email or archive.")),
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
        card_header(bsicons::bs_icon("layers"),
                    tags$strong(" MacTel gene tiers"),
                    tags$span(" — the curated candidate-gene list",
                              class = "text-body-secondary")),
        tags$p(class = "text-body-secondary small mb-3",
          "Genes are grouped by how strong the evidence is that they are ",
          "involved in MacTel. Use the ", tags$strong("Gene Tier"),
          " filter on the left to focus on one tier at a time."),
        layout_columns(
          col_widths = c(6, 6),
          tier_box("Tier 1", textOutput("n_tier1", inline = TRUE),
                   "Established, high-confidence links to MacTel.",
                   list(bg = "#eaf1f8", bar = "#1F4E79")),
          tier_box("Tier 2", textOutput("n_tier2", inline = TRUE),
                   "Candidate genes with supporting but less definitive evidence.",
                   list(bg = "#eef5f9", bar = "#3A7CA5"))
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
      div(class = "d-flex justify-content-end mb-2",
          actionButton("open_dl_overview", "Download figure (PNG)",
                       class = "btn-sm btn-outline-secondary",
                       icon = bsicons::bs_icon("image"))),
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
                         class = "btn-sm btn-primary float-end"),
          actionButton("xl_table", "Export to Excel",
                       icon = bsicons::bs_icon("file-earmark-spreadsheet"),
                       class = "btn-sm btn-success float-end me-2")
        ),
        uiOutput("legend_variants"),
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
                         class = "btn-sm btn-primary float-end"),
          actionButton("xl_priority", "Export to Excel",
                       icon = bsicons::bs_icon("file-earmark-spreadsheet"),
                       class = "btn-sm btn-success float-end me-2")
        ),
        uiOutput("legend_priority"),
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
      div(
        # The sample picker sits in a compact header strip directly under the
        # tab bar (rather than a left sidebar), so the content below can use the
        # full width instead of leaving a tall empty column.
        div(
          class = "mb-3",
          style = "min-width: 280px; max-width: 360px;",
          selectizeInput("sample_pick", "Select a sample", width = "100%",
                         choices = NULL, multiple = FALSE,
                         options = list(placeholder = "Start typing a sample ID…"))
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
          ),
          card(
            card_header(
              "IGV report",
              tags$span(bsicons::bs_icon("info-circle"),
                        " read-level view of every variant this sample carries ",
                        "(requires an internet connection)",
                        class = "text-muted small ms-2")
            ),
            uiOutput("sample_igv")
          )
        )
      )
    ),

    nav_panel(
      "Family explorer",
      icon = bsicons::bs_icon("people-fill"),
      div(
        div(
          class = "mb-3",
          style = "min-width: 280px; max-width: 360px;",
          selectizeInput("family_pick", "Select a family", width = "100%",
                         choices = c("", FAMILY_CHOICES), selected = "",
                         options = list(placeholder = "Choose a family…"))
        ),
        uiOutput("family_header"),
        uiOutput("family_body")
      )
    )
  )
)

# -----------------------------------------------------------------------------
# SERVER
# -----------------------------------------------------------------------------
# Tracks how many browser sessions are connected, so the desktop launcher can
# quit the R process (and free the port) once the last tab is closed.
.autostop <- new.env()
.autostop$n <- 0L

server <- function(input, output, session) {

  # When launched via the desktop launcher (which sets this option), shut the
  # app down shortly after the last browser tab closes — this frees port 7766
  # so the next launch starts cleanly. A 2-second grace period means a page
  # refresh (old session ends, new one connects) does NOT trigger a shutdown.
  if (isTRUE(getOption("mactel.autostop", FALSE))) {
    .autostop$n <- .autostop$n + 1L
    session$onSessionEnded(function() {
      .autostop$n <- .autostop$n - 1L
      quit_if_idle <- function() if (.autostop$n <= 0L) stopApp()
      if (requireNamespace("later", quietly = TRUE))
        later::later(quit_if_idle, delay = 2)
      else
        quit_if_idle()
    })
  }

  # ---- raw data (reactive on upload) ----------------------------------------
  raw <- reactiveVal(NULL)
  src_label <- reactiveVal("No file loaded — upload a Cavalier CSV to begin")
  main_name <- reactiveVal("")                       # label for the active file

  # ---- active gene-list version --------------------------------------------
  # The selected version supplies both the gene descriptions (gene_info_rv) and
  # the SYMBOL -> Tier lookup (tier_df_rv). Switching versions re-tiers whatever
  # variants are loaded and refreshes the tier filter choices downstream.
  gene_info_rv <- reactiveVal(GENE_INFO)
  tier_df_rv   <- reactiveVal(TIER_DF)

  observeEvent(input$gene_list_pick, {
    path <- GENE_LISTS[[input$gene_list_pick]]
    if (is.null(path) || !file.exists(path)) return()
    gi  <- load_gene_info(path)
    tdf <- if (!is.null(gi))
      dplyr::distinct(dplyr::select(gi, SYMBOL, Tier)) else NULL
    gene_info_rv(gi)
    tier_df_rv(tdf)
    # Re-tier the currently loaded variants under the new gene list.
    df <- raw()
    if (!is.null(df)) {
      df <- dplyr::select(df, -dplyr::any_of("Tier"))
      raw(annotate_tier(df, tdf))
    }
  }, ignoreInit = TRUE)

  # Landing-page tier counts track the active gene list.
  output$n_tier1 <- renderText(tier_gene_count(tier_df_rv(), "Tier 1"))
  output$n_tier2 <- renderText(tier_gene_count(tier_df_rv(), "Tier 2"))

  # Auto-load the startup file only in DEBUG mode (developer convenience). Normal
  # launches stay empty so the user explicitly chooses their own file.
  if (isTRUE(DEBUG) && file.exists(startup_path)) {
    main_name(basename(startup_path))
    observe({
      df <- load_annotated(startup_path, isolate(tier_df_rv()))
      raw(df)
      src_label(sprintf("%s  (%d variants, DEBUG auto-load)",
                        basename(startup_path), nrow(df)))
    })
  }

  observeEvent(input$upload, {
    req(input$upload)
    df <- tryCatch(load_annotated(input$upload$datapath, tier_df_rv()),
                   error = function(e) {
                     showNotification(paste("Could not load file:", e$message),
                                      type = "error", duration = 8)
                     NULL
                   })
    if (!is.null(df)) {
      raw(df)
      main_name(input$upload$name)
      src_label(sprintf("%s  (%d variants, uploaded)",
                        input$upload$name, nrow(df)))
    }
  })

  output$data_source_label <- renderText(paste("Source:", src_label()))

  # ---- variant comparison (second uploaded file) ----------------------------
  # Uploading a second CSV opens a "Variant comparisons" tab listing variants
  # present in only one of the two files (compared on CHROM:POS:REF:ALT).
  compare_raw       <- reactiveVal(NULL)
  compare_name      <- reactiveVal(NULL)
  compare_tab_added <- reactiveVal(FALSE)

  # The family whose members are shown in the Family explorer tab. NULL until a
  # clickable Family_ID badge is used in the Sample explorer.
  selected_family <- reactiveVal(NULL)

  observeEvent(input$compare_upload, {
    req(input$compare_upload)
    df <- tryCatch(load_annotated(input$compare_upload$datapath, tier_df_rv()),
                   error = function(e) {
                     showNotification(paste("Could not load comparison file:",
                                            e$message),
                                      type = "error", duration = 8)
                     NULL
                   })
    req(!is.null(df))
    compare_raw(df)
    compare_name(input$compare_upload$name)

    if (!isTRUE(compare_tab_added())) {
      bslib::nav_insert(
        id = "main_tabs",
        target = "Variant table", position = "after", select = TRUE,
        session = session,
        nav = nav_panel(
          "Variant comparisons",
          icon = bsicons::bs_icon("intersect"),
          card(
            card_header(
              "Variants unique to one file",
              tags$span(bsicons::bs_icon("info-circle"),
                        textOutput("compare_summary", inline = TRUE),
                        class = "text-muted small ms-2"),
              downloadButton("dl_compare", "Download CSV",
                             class = "btn-sm btn-primary float-end")
            ),
            DT::DTOutput("compare_table")
          )
        )
      )
      compare_tab_added(TRUE)
    } else {
      bslib::nav_select("main_tabs", "Variant comparisons", session = session)
    }
  })

  # Variants in exactly one file, de-duplicated to one row per variant, tagged
  # with which file they are unique to.
  compare_diff <- reactive({
    a <- raw(); b <- compare_raw()
    req(a, b)
    ak <- paste(a$CHROM, a$POS, a$REF, a$ALT)
    bk <- paste(b$CHROM, b$POS, b$REF, b$ALT)
    a_only <- a[!(ak %in% bk), , drop = FALSE] %>%
      dplyr::distinct(CHROM, POS, REF, ALT, .keep_all = TRUE) %>%
      dplyr::mutate(cmp_src = main_name())
    b_only <- b[!(bk %in% ak), , drop = FALSE] %>%
      dplyr::distinct(CHROM, POS, REF, ALT, .keep_all = TRUE) %>%
      dplyr::mutate(cmp_src = compare_name())
    dplyr::bind_rows(a_only, b_only)
  })

  display_compare <- function(df, links = TRUE) {
    gene_col    <- if (links) link_gene(df$SYMBOL) else df$SYMBOL
    variant_col <- if (links) link_variant(df$CHROM, df$POS, df$REF, df$ALT)
                   else sprintf("%s:%s %s>%s", df$CHROM, df$POS, df$REF, df$ALT)
    df %>%
      dplyr::transmute(
        `Only in` = cmp_src,
        Gene = gene_col, Tier = Tier,
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

  output$compare_summary <- renderText({
    d <- compare_diff(); req(d)
    na <- sum(d$cmp_src == main_name())
    nb <- sum(d$cmp_src == compare_name())
    sprintf(" %d only in %s · %d only in %s",
            na, main_name(), nb, compare_name())
  })

  output$compare_table <- DT::renderDT({
    d <- compare_diff()
    validate(need(nrow(d) > 0,
                  "No differences — both files contain the same variants."))
    dt <- display_compare(d)
    DT::datatable(dt,
                  filter = "top", rownames = FALSE,
                  selection = "none", escape = FALSE,
                  extensions = "Buttons",
                  options = list(pageLength = 25, scrollX = TRUE,
                                 order = order_desc_by(dt, "CADD"),
                                 dom = "Bfrtip", buttons = c("copy", "csv"),
                                 headerCallback = header_tips_cb())) %>%
      DT::formatStyle("Only in", fontWeight = "bold") %>%
      DT::formatStyle("ClinVar",
                      backgroundColor = DT::styleEqual(
                        c("Pathogenic", "Pathogenic/Likely_pathogenic",
                          "Likely_pathogenic"),
                        c("#FDDEDE", "#F7D6F7", "#FDE8D8"))) %>%
      DT::formatStyle("Impact",
                      backgroundColor = DT::styleEqual("HIGH", "#FDDEDE"))
  })

  output$dl_compare <- downloadHandler(
    filename = function() sprintf("MacTel_variant_comparison_%s.csv", Sys.Date()),
    content  = function(file)
      readr::write_csv(display_compare(compare_diff(), links = FALSE), file)
  )

  # ---- populate dynamic filter choices when data changes --------------------
  observeEvent(raw(), {
    df <- raw()
    tiers <- sort(unique(df$Tier))
    updateCheckboxGroupInput(session, "tier",
                             choices = tiers, selected = tiers, inline = TRUE)
    updateSelectizeInput(session, "genes",
                         choices = sort(unique(df$SYMBOL)), server = TRUE)
    fids <- sort(unique(df$family_id))
    updateSelectizeInput(session, "sample_pick",
                         choices = stats::setNames(fids, fmt_sample(fids)),
                         server = TRUE)
    updateSelectizeInput(session, "exclude_samples",
                         choices = stats::setNames(fids, fmt_sample(fids)),
                         selected = isolate(input$exclude_samples) %||% character(0),
                         server = TRUE)
    mx <- ceiling(max(df$CADD, na.rm = TRUE))
    updateSliderInput(session, "cadd", max = mx, value = 0)
    updateSliderInput(session, "priority_cadd", max = mx)
  })

  # Re-label the sample picker when the ID format toggle changes (the value
  # stays family_id; only the visible label changes). Tables/plots re-render
  # on their own because they read input$id_format via fmt_sample().
  observeEvent(input$id_format, {
    df <- raw(); req(df)
    fids <- sort(unique(df$family_id))
    updateSelectizeInput(session, "sample_pick",
                         choices = stats::setNames(fids, fmt_sample(fids)),
                         selected = isolate(input$sample_pick) %||% "",
                         server = TRUE)
    updateSelectizeInput(session, "exclude_samples",
                         choices = stats::setNames(fids, fmt_sample(fids)),
                         selected = isolate(input$exclude_samples) %||% character(0),
                         server = TRUE)
  }, ignoreInit = TRUE)

  observeEvent(input$reset_filters, {
    df <- raw(); req(df)
    updateCheckboxGroupInput(session, "tier",
                             selected = sort(unique(df$Tier)))
    updateSelectizeInput(session, "genes", selected = character(0))
    updateCheckboxGroupInput(session, "impact", selected = IMPACT_LEVELS)
    updateCheckboxGroupInput(session, "type", selected = TYPE_LEVELS)
    updateCheckboxGroupInput(session, "clnsig", selected = CLNSIG_LEVELS)
    updateCheckboxInput(session, "exclude_am_benign", value = FALSE)
    updateSliderInput(session, "cadd", value = 0)
    updateSliderInput(session, "revel", value = 0)
    updateNumericInput(session, "gnomad", value = NA)
    updateSliderInput(session, "min_flags", value = 0)
    updateCheckboxGroupInput(session, "sample_group", selected = "MacTel")
    updateSelectizeInput(session, "exclude_samples", selected = character(0))
  })

  # ---- core filtered dataset ------------------------------------------------
  # All filters EXCEPT the sample-group filter. The sample explorer uses this so
  # an individual can be searched even if their diagnosis group is unticked.
  filtered_pre_group <- reactive({
    df <- raw(); req(df)

    # Global sample exclusion — drops chosen samples from every view,
    # independent of the sample-group selection below.
    if (length(input$exclude_samples) > 0)
      df <- dplyr::filter(df, !family_id %in% input$exclude_samples)

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
    if (isTRUE(input$exclude_am_benign))
      df <- dplyr::filter(df, is.na(am_class) |
                              !tolower(am_class) %in% c("benign", "likely_benign"))

    df <- dplyr::filter(df, is.na(CADD) | CADD >= input$cadd)
    if (input$revel > 0)
      df <- dplyr::filter(df, !is.na(REVEL) & REVEL >= input$revel)
    if (!is.null(input$gnomad) && !is.na(input$gnomad) && input$gnomad < 1) {
      thr <- input$gnomad
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
  output$p_genes   <- renderPlot(
    plot_top_genes(filtered(), 25, group_lookup = DIAG_GROUP_LOOKUP))

  # ---- overview figure export (multi-panel PNG of the current filters) ------
  # Concise one-line description of which filters are currently active, so the
  # exported figure is self-documenting. Only non-default filters are listed.
  overview_caption <- reactive({
    d <- filtered()
    counts <- sprintf(
      "%s variants | %d genes | %d samples | %d ClinVar P/LP",
      format(nrow(d), big.mark = ","), dplyr::n_distinct(d$SYMBOL),
      dplyr::n_distinct(d$family_id), sum(d$is_pathLP))
    parts <- c()
    all_tiers <- sort(unique(raw()$Tier))
    if (length(input$tier) && !setequal(input$tier, all_tiers))
      parts <- c(parts, paste("Tier", paste(sort(input$tier), collapse = ",")))
    if (length(input$genes))
      parts <- c(parts, paste("Genes", paste(input$genes, collapse = ",")))
    if (length(input$impact) && !setequal(input$impact, IMPACT_LEVELS))
      parts <- c(parts, paste("Impact", paste(input$impact, collapse = ",")))
    if (length(input$type) && !setequal(input$type, TYPE_LEVELS))
      parts <- c(parts, paste("Type", paste(input$type, collapse = ",")))
    if (length(input$clnsig) && !setequal(input$clnsig, CLNSIG_LEVELS))
      parts <- c(parts, paste("ClinVar", paste(input$clnsig, collapse = ",")))
    if (!is.null(input$cadd) && input$cadd > 0)
      parts <- c(parts, sprintf("CADD>=%g", input$cadd))
    if (!is.null(input$revel) && input$revel > 0)
      parts <- c(parts, sprintf("REVEL>=%g", input$revel))
    if (!is.null(input$gnomad) && !is.na(input$gnomad) && input$gnomad < 1)
      parts <- c(parts, sprintf("gnomAD<=%g", input$gnomad))
    if (isTRUE(input$exclude_am_benign))
      parts <- c(parts, "AlphaMissense benign excluded")
    if (!is.null(SAMPLE_INFO) && length(input$sample_group) < 3)
      parts <- c(parts, paste("Group", paste(input$sample_group, collapse = ",")))
    filt <- if (length(parts)) paste(parts, collapse = " | ") else
      "no filters applied (all variants)"
    paste0(counts, "\nFilters: ", filt)
  })

  # Panel keys -> human labels for the download picker. Order here is the
  # order panels are laid out in the exported figure.
  OVERVIEW_PANELS <- c(
    impact  = "VEP impact",
    type    = "Variant type",
    inherit = "Inheritance mode",
    cadd    = "CADD distribution",
    clnsig  = "ClinVar classification",
    genes   = "Top genes by samples"
  )

  # Colour-palette options for the figure export. Labels are shown to the user;
  # values are passed to apply_palette() ("Default"/"Colour-blind" are handled
  # specially, everything else is a grDevices::hcl.colors palette name).
  OVERVIEW_PALETTES <- c(
    "Default (semantic)"        = "Default",
    "Viridis"                   = "Viridis",
    "Cividis (colour-blind)"    = "Cividis",
    "Plasma"                    = "Plasma",
    "Inferno"                   = "Inferno",
    "Mako"                      = "Mako",
    "Okabe-Ito (colour-blind)"  = "Colour-blind",
    "Set 2 (soft)"              = "Set 2",
    "Dark 2 (bold)"             = "Dark 2",
    "Dark 3"                    = "Dark 3",
    "Zissou 1"                  = "Zissou 1",
    "Temps"                     = "Temps"
  )

  # Ask the user which panels + colour palette before generating the figure.
  observeEvent(input$open_dl_overview, {
    showModal(modalDialog(
      title = "Download overview figure",
      size  = "m",
      easyClose = TRUE,
      checkboxGroupInput(
        "ov_panels", "Panels to include:",
        choices  = stats::setNames(names(OVERVIEW_PANELS), unname(OVERVIEW_PANELS)),
        selected = names(OVERVIEW_PANELS)),
      selectInput(
        "ov_palette", "Colour palette:",
        choices  = OVERVIEW_PALETTES,
        selected = "Default"),
      footer = tagList(
        modalButton("Cancel"),
        downloadButton("dl_overview", "Download PNG",
                       class = "btn-primary",
                       icon = bsicons::bs_icon("download")))
    ))
  })

  output$dl_overview <- downloadHandler(
    filename = function() sprintf("mactel_overview_%s.png", Sys.Date()),
    content  = function(file) {
      d   <- filtered()
      pal <- input$ov_palette %||% "Default"
      sel <- input$ov_panels
      if (is.null(sel) || length(sel) == 0) sel <- names(OVERVIEW_PANELS)
      # keep the canonical panel order regardless of tick order
      sel <- names(OVERVIEW_PANELS)[names(OVERVIEW_PANELS) %in% sel]

      # A blank placeholder keeps the grid intact when a plot builder returns
      # NULL (e.g. no CADD values, or no genes after filtering).
      or_blank <- function(p, msg) if (!is.null(p)) p else
        ggplot2::ggplot() +
          ggplot2::annotate("text", x = 0, y = 0, label = msg,
                            colour = "grey50", size = 4) +
          ggplot2::theme_void()

      builders <- list(
        impact  = function() plot_impact(d, palette = pal),
        type    = function() plot_type(d, palette = pal),
        inherit = function() plot_inheritance(d, palette = pal),
        cadd    = function() or_blank(
          plot_cadd(d, input$priority_cadd, palette = pal), "No CADD values"),
        clnsig  = function() plot_clnsig(d, palette = pal),
        genes   = function() or_blank(
          plot_top_genes(d, 15, group_lookup = DIAG_GROUP_LOOKUP, palette = pal),
          "No genes to show")
      )
      panels <- lapply(sel, function(k) builders[[k]]())

      # "genes" is a tall horizontal bar chart, so give it a full-width row of
      # its own; the remaining panels flow two-per-row above it.
      has_genes  <- "genes" %in% sel
      compact    <- panels[sel != "genes"]
      comp_rows  <- if (length(compact) > 0) ceiling(length(compact) / 2) else 0
      # Physical height (inches) for each block, so the compact grid keeps a
      # constant per-row height instead of being squashed into one slot.
      COMPACT_ROW_H <- 3.3
      GENES_H       <- 5.5
      compact_h  <- comp_rows * COMPACT_ROW_H
      genes_h    <- if (has_genes) GENES_H else 0

      pieces  <- list()
      heights <- numeric(0)
      if (length(compact) > 0) {
        pieces  <- c(pieces, list(patchwork::wrap_plots(compact, ncol = 2)))
        heights <- c(heights, compact_h)
      }
      if (has_genes) {
        pieces  <- c(pieces, list(panels[[which(sel == "genes")]]))
        heights <- c(heights, genes_h)
      }

      fig <- patchwork::wrap_plots(pieces, ncol = 1, heights = heights) +
        patchwork::plot_annotation(
          title    = "MacTel Variant Explorer - overview",
          subtitle = overview_caption(),
          caption  = format(Sys.Date()),
          theme = ggplot2::theme(
            plot.title    = ggplot2::element_text(face = "bold", size = 18),
            plot.subtitle = ggplot2::element_text(size = 11, colour = "grey30"),
            plot.caption  = ggplot2::element_text(size = 9, colour = "grey50")))

      # Total canvas height matches the sum of block heights (+ a little for the
      # title/subtitle), so nothing is stretched or crushed.
      fig_h <- max(5, compact_h + genes_h + 1)
      # device = "png" is required: downloadHandler hands us an extension-less
      # temp path, so ggsave cannot infer the format from the filename.
      ggplot2::ggsave(file, fig, device = "png", width = 13,
                      height = fig_h, dpi = 200, bg = "white")
      removeModal()
    }
  )

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
  # Display label for a family_id under the currently selected ID format.
  # Reading input$id_format here makes every table/plot that calls it
  # re-render when the user flips the format toggle.
  fmt_sample <- function(fid) format_sample_id(fid, input$id_format %||% "AID")
  link_gene <- function(sym) {
    ifelse(is.na(sym) | sym == "", as.character(sym), sprintf(
      "<a href='#' onclick=\"Shiny.setInputValue('cell_gene','%s',{priority:'event'});return false;\">%s</a>",
      .jsesc(sym), sym))
  }
  # Payload stays family_id (so jump-to-sample keeps working); the visible text
  # uses the selected ID format. Rendered as a diagnosis-coloured pill badge so
  # a sample's MacTel/HSAN1/both/control status is readable straight from the
  # table without having to click through. Vectorised over fid.
  # onclick (carrying family_id) is placed before the diagnosis-dependent style
  # so the DT sort key prefix stays constant and the column still sorts by id.
  link_sample <- function(fid) sprintf(
    "<a href='#' onclick=\"Shiny.setInputValue('cell_sample','%s',{priority:'event'});return false;\" class='badge rounded-pill me-1' style='background-color:%s;color:#fff;text-decoration:none;cursor:pointer;'>%s</a>",
    .jsesc(fid), diag_colour(fid), fmt_sample(fid))
  # input_id lets the in-modal variant table fire a different input
  # (modal_pick_variant) so it toggles state in place instead of re-opening.
  link_variant <- function(chrom, pos, ref, alt, input_id = "cell_variant") {
    label <- sprintf("%s:%s %s>%s", chrom, pos, ref, alt)
    key   <- sprintf("%s||%s||%s||%s", chrom, pos, ref, alt)
    sprintf(
      "<a href='#' onclick=\"Shiny.setInputValue('%s','%s',{priority:'event'});return false;\">%s</a>",
      input_id, .jsesc(key), label)
  }
  # Diagnosis colour for one or more samples (family_id keys): red MacTel, blue
  # HSAN1, purple both, grey control/unknown. Vectorised (via DIAG_GROUP_LOOKUP +
  # COL_DIAG) so it can colour a whole table column at once. Shared by the sample
  # badges, the variant-modal carrier chips, and the Family explorer.
  diag_colour <- function(fid) {
    fid <- as.character(fid)
    if (is.null(DIAG_GROUP_LOOKUP)) return(rep("#546E7A", length(fid)))
    col <- unname(COL_DIAG[unname(DIAG_GROUP_LOOKUP[fid])])
    col[is.na(col)] <- "#546E7A"
    col
  }
  # Diagnosis group label(s) for one or more samples (for legends/tooltips).
  diag_group <- function(fid) {
    fid <- as.character(fid)
    if (is.null(DIAG_GROUP_LOOKUP)) return(rep(NA_character_, length(fid)))
    unname(DIAG_GROUP_LOOKUP[fid])
  }
  # link_member is kept as an alias so existing callers still work; sample links
  # everywhere now render as the same clickable, diagnosis-coloured badge.
  link_member <- function(fid) link_sample(fid)

  # Small inline legend of the diagnosis colours actually present among `fids`,
  # shown above sample-bearing tables so the badge colours are self-explanatory.
  # Returns NULL when there is no sample info or no group to show.
  diag_legend <- function(fids) {
    if (is.null(DIAG_GROUP_LOOKUP)) return(NULL)
    present <- intersect(names(COL_DIAG),
                         unique(diag_group(fids)[!is.na(diag_group(fids))]))
    if (length(present) == 0) return(NULL)
    swatch <- function(lab) tags$span(
      class = "badge rounded-pill me-2",
      style = sprintf("background-color:%s;color:#fff;font-weight:normal;",
                      COL_DIAG[[lab]]), lab)
    tags$div(class = "mb-2 small",
             tags$span("Sample diagnosis:", class = "text-muted me-2"),
             lapply(present, swatch))
  }

  # DT `order` spec that sorts a table by the first matching column name,
  # descending, by default. Column index is 0-based (rownames are off in every
  # table here). Falls back to DT's default ordering when no column matches.
  order_desc_by <- function(df, cols) {
    idx <- match(intersect(cols, names(df))[1], names(df))
    if (is.na(idx)) list() else list(list(idx - 1L, "desc"))
  }

  # ---- display-table builder ------------------------------------------------
  display_cols <- function(df, links = TRUE) {
    if (links) {
      gene_col    <- link_gene(df$SYMBOL)
      sample_col  <- link_sample(df$family_id)
      variant_col <- link_variant(df$CHROM, df$POS, df$REF, df$ALT)
    } else {
      gene_col    <- df$SYMBOL
      sample_col  <- fmt_sample(df$family_id)
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

  output$legend_variants <- renderUI(diag_legend(filtered()$family_id))

  output$variant_table <- DT::renderDT({
    dt <- display_cols(filtered())
    DT::datatable(dt,
                  filter = "top", rownames = FALSE,
                  selection = "none", escape = FALSE,
                  extensions = "Buttons",
                  options = list(pageLength = 25, scrollX = TRUE,
                                 order = order_desc_by(dt, "CADD"),
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
  output$legend_priority <- renderUI(diag_legend(priority()$family_id))

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
                                 order = order_desc_by(tbl, "CADD"),
                                 headerCallback = header_tips_cb())) %>%
      DT::formatStyle("Flags", fontWeight = "bold",
                      background = DT::styleColorBar(c(0, 3), "#9ec5fe")) %>%
      DT::formatStyle("Impact",
                      backgroundColor = DT::styleEqual("HIGH", "#FDDEDE"))
  })

  output$p_priority_genes <- renderPlot({
    d <- priority()
    validate(need(nrow(d) > 0, "No variants meet the chosen number of flags."))
    plot_top_genes(d, 20, group_lookup = DIAG_GROUP_LOOKUP)
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
      dplyr::arrange(dplyr::desc(Variants), dplyr::desc(Samples))
  })

  output$gene_table <- DT::renderDT({
    g <- gene_summary() %>% dplyr::mutate(SYMBOL = link_gene(SYMBOL))
    DT::datatable(g, rownames = FALSE, filter = "top",
                  selection = "none", escape = FALSE,
                  options = list(pageLength = 25, scrollX = TRUE,
                                 order = order_desc_by(g, "Variants"),
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
    else sprintf("Variants in sample %s", fmt_sample(input$sample_pick))
  })

  # ---- family explorer ------------------------------------------------------
  # All sequenced members of the selected family (family_id keys), used as the
  # denominator for carrier counts and the universe for both family tables.
  family_member_keys <- reactive({
    fam <- selected_family()
    if (is.null(fam) || is.null(FAMILY_MEMBERS[[fam]])) return(character(0))
    FAMILY_MEMBERS[[fam]]
  })

  family_data_priority <- reactive({
    members <- family_member_keys(); req(length(members) > 0)
    filtered_pre_group() %>% dplyr::filter(family_id %in% members)
  })

  family_data_all <- reactive({
    members <- family_member_keys(); req(length(members) > 0)
    raw() %>% dplyr::filter(family_id %in% members)
  })

  # Collapse a family's variant rows to one row per unique variant, adding a
  # carrier count (N members carrying / total sequenced) and a clickable,
  # diagnosis-coloured member-badge column. Sorted by carrier count desc.
  build_family_dt <- function(d, total) {
    if (is.null(d) || nrow(d) == 0) return(NULL)
    d <- d %>% dplyr::mutate(
      .vkey = paste(CHROM, POS, REF, ALT, sep = "||"))
    carriers <- d %>%
      dplyr::distinct(.vkey, family_id) %>%
      dplyr::group_by(.vkey) %>%
      dplyr::summarise(
        n_carriers = dplyr::n_distinct(family_id),
        Members = paste(vapply(sort(unique(family_id)), link_member,
                               character(1)), collapse = " "),
        .groups = "drop")
    reps <- d %>%
      dplyr::group_by(.vkey) %>% dplyr::slice(1) %>% dplyr::ungroup()
    tbl <- reps %>%
      dplyr::left_join(carriers, by = ".vkey") %>%
      dplyr::arrange(dplyr::desc(n_carriers), dplyr::desc(is_pathLP),
                     dplyr::desc(CADD)) %>%
      dplyr::transmute(
        Gene = link_gene(SYMBOL), Tier = Tier,
        Variant = link_variant(CHROM, POS, REF, ALT),
        Carriers = sprintf("%d/%d", n_carriers, total),
        Members = Members,
        HGVSc, HGVSp = HGVSp_short,
        Impact = IMPACT, Type = TYPE, CADD = round(CADD, 1),
        REVEL = round(REVEL, 3), AlphaMissense = am_class,
        ClinVar = CLNSIG_clean, gnomAD_AF = signif(gnomad_AF, 3),
        Inheritance = inheritance)
    DT::datatable(tbl, rownames = FALSE, selection = "none", escape = FALSE,
                  # Rows arrive pre-sorted by descending carrier count; an empty
                  # order preserves that (the "n/total" Carriers string cannot be
                  # sorted numerically by DT).
                  options = list(pageLength = 15, scrollX = TRUE,
                                 order = list(),
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

  output$family_header <- renderUI({
    fam <- selected_family()
    if (is.null(fam)) return(NULL)
    members <- family_member_keys()
    legend_dot <- function(col, lab) tags$span(
      tags$span(style = sprintf(
        "display:inline-block;width:10px;height:10px;border-radius:50%%;background-color:%s;margin-right:3px;",
        col)), lab, class = "me-2")
    tagList(
      div(class = "mb-2 d-flex flex-wrap align-items-center",
          tags$span(class = "me-3",
            tags$span("Family ", class = "text-muted small"),
            tags$span(fam, class = "fw-semibold")),
          tags$span(
            tags$span("Sequenced members ", class = "text-muted small"),
            tags$span(length(members), class = "fw-semibold"))),
      div(class = "mb-2", HTML(paste(vapply(members, link_member,
                                            character(1)), collapse = " "))),
      div(class = "small text-muted mb-3",
          legend_dot("#C62828", "MacTel"), legend_dot("#1565C0", "HSAN1"),
          legend_dot("#6A1B9A", "MacTel + HSAN1"),
          legend_dot("#546E7A", "Control"))
    )
  })

  output$family_body <- renderUI({
    if (is.null(selected_family()))
      return(div(class = "text-muted",
                 bsicons::bs_icon("people"),
                 paste0(" Choose a family above, or click a Family badge in the ",
                        "Sample explorer, to see the whole family's variants ",
                        "here.")))
    tagList(
      card(
        card_header(
          "Filtered / prioritised variants",
          tags$span(bsicons::bs_icon("info-circle"),
                    " family variants that pass the global filters; one row per ",
                    "unique variant",
                    class = "text-muted small ms-2")
        ),
        DT::DTOutput("family_table_priority")
      ),
      card(
        card_header(
          "All variants",
          tags$span(bsicons::bs_icon("info-circle"),
                    " every variant carried by any family member, ignoring ",
                    "filters; one row per unique variant",
                    class = "text-muted small ms-2")
        ),
        DT::DTOutput("family_table_all")
      )
    )
  })

  output$family_table_priority <- DT::renderDT({
    validate(need(!is.null(selected_family()),
                  "Click a Family badge in the Sample explorer."))
    d <- family_data_priority()
    validate(need(nrow(d) > 0,
                  "No family variants pass the current filters."))
    build_family_dt(d, length(family_member_keys()))
  })

  output$family_table_all <- DT::renderDT({
    validate(need(!is.null(selected_family()),
                  "Click a Family badge in the Sample explorer."))
    d <- family_data_all()
    validate(need(nrow(d) > 0, "No variants for this family."))
    build_family_dt(d, length(family_member_keys()))
  })

  # ---- per-sample IGV report ------------------------------------------------
  # Point the app at a folder of igv-reports HTMLs (one sub-folder per sample)
  # using a point-and-click folder picker (shinyFiles). The chosen folder is
  # exposed over HTTP via a resource path so it can be shown in an <iframe>.
  # The display family_id has its trailing "RLA" stripped, but the folders keep
  # it, so we probe both the displayed ID and "<id>RLA".
  igv_volumes <- c(Home = path.expand("~"), shinyFiles::getVolumes()())

  # Register a folder and remember it; returns the normalised path (or NULL).
  use_igv_root <- function(d) {
    if (length(d) == 1 && nzchar(d) && dir.exists(d)) {
      root <- normalizePath(d)
      suppressWarnings(shiny::addResourcePath("igv_reports", root))
      igv_root(root)
      return(root)
    }
    NULL
  }

  # Seed with the auto-detected folder (data/igv_reports) when present.
  igv_root <- reactiveVal(NULL)
  isolate(if (nzchar(igv_startup) && dir.exists(igv_startup))
            use_igv_root(igv_startup))

  shinyFiles::shinyDirChoose(input, "igv_dir_btn", roots = igv_volumes,
                             session = session, allowDirCreate = FALSE)
  observeEvent(input$igv_dir_btn, {
    req(is.list(input$igv_dir_btn))
    d <- shinyFiles::parseDirPath(igv_volumes, input$igv_dir_btn)
    use_igv_root(d)
  })

  output$igv_dir_label <- renderText({
    root <- igv_root()
    if (is.null(root)) "No folder selected" else root
  })

  # Resolve the on-disk folder name for the picked sample (NULL if no report).
  igv_sample_dir <- function(root, sample) {
    for (s in unique(c(sample, paste0(sample, "RLA")))) {
      f <- file.path(root, s, paste0(s, ".igv_report.html"))
      if (file.exists(f)) return(s)
    }
    NULL
  }

  output$sample_igv <- renderUI({
    note <- function(msg)
      tags$p(class = "text-muted small mb-0", msg)

    if (is.null(input$sample_pick) || input$sample_pick == "")
      return(note("Choose a sample from the dropdown."))
    root <- igv_root()
    if (is.null(root))
      return(note(paste("Set the 'IGV reports folder' in the Data panel",
                        "to view read-level IGV reports here.")))

    hit <- igv_sample_dir(root, input$sample_pick)
    if (is.null(hit))
      return(note(sprintf("No IGV report found for sample %s.",
                          input$sample_pick)))

    src <- sprintf("igv_reports/%s/%s.igv_report.html",
                   utils::URLencode(hit, reserved = TRUE),
                   utils::URLencode(hit, reserved = TRUE))
    tags$iframe(
      src = src, loading = "lazy",
      style = paste("width:100%; height:680px; border:1px solid #dee2e6;",
                    "border-radius:4px;"))
  })

  # Data-group tags for the picked sample (Mito_haplo deliberately excluded).
  output$sample_tags <- renderUI({
    if (is.null(SAMPLE_INFO) ||
        is.null(input$sample_pick) || input$sample_pick == "")
      return(NULL)
    row <- SAMPLE_INFO[SAMPLE_INFO$family_id == input$sample_pick, , drop = FALSE]
    if (nrow(row) == 0) return(NULL)
    row <- row[1, ]

    # Diagnosis badge first. Colours match the carrier/member badges elsewhere
    # (MacTel red, HSAN1 blue, control slate) with white text for contrast.
    diag_pill <- function(label, bg) tags$span(
      label, class = "badge rounded-pill me-1",
      style = sprintf("background-color:%s;color:#fff;", bg))
    diag_badge <- if (isTRUE(row$is_mactel)) {
      diag_pill("MacTel", "#C62828")
    } else if (isTRUE(row$is_hsan1)) {
      diag_pill("HSAN1", "#1565C0")
    } else {
      diag_pill("Control", "#546E7A")
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

    # Identity line: always show BOTH identifiers regardless of the selected
    # display format. Prefer the non-padded AID form (e.g. A1) from the sample
    # sheet; fall back to the family_id key if the column is missing.
    aid <- if ("AID" %in% names(row) &&
               !is.na(row$AID) && nzchar(row$AID))
      row$AID else input$sample_pick
    pid <- if ("Patient_ID" %in% names(row) &&
               !is.na(row$Patient_ID) && nzchar(row$Patient_ID))
      row$Patient_ID else NA_character_
    fam <- family_of(input$sample_pick)
    fam_badge <- if (!is.na(fam)) {
      tags$span(class = "ms-3",
        tags$span("Family ", class = "text-muted small"),
        tags$a(href = "#", fam,
          class = "fw-semibold badge rounded-pill",
          style = paste0("background-color:#00695C;color:#fff;",
                         "text-decoration:none;cursor:pointer;"),
          onclick = sprintf(
            "Shiny.setInputValue('cell_family','%s',{priority:'event'});return false;",
            .jsesc(fam))))
    }
    id_line <- div(
      class = "mb-2 d-flex flex-wrap align-items-center",
      tags$span(class = "me-3",
        tags$span("AID ", class = "text-muted small"),
        tags$span(aid, class = "fw-semibold")),
      if (!is.na(pid)) tags$span(
        tags$span("Patient ID ", class = "text-muted small"),
        tags$span(pid, class = "fw-semibold")),
      fam_badge
    )

    tagList(
      id_line,
      div(class = "mb-2",
          tags$span("Tags: ", class = "text-muted small me-1"),
          diag_badge, group_badges)
    )
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
                                 order = order_desc_by(tbl, "CADD"),
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

  # ---- portable filter codes (share / save / restore) -----------------------
  # Instead of a URL that bundles the whole app state (and needs the same data
  # file to make sense), we serialise just the sidebar filter settings into a
  # short base64 code. A colleague pastes that code to apply the identical
  # filters to their OWN loaded data, and the same code doubles as a "save my
  # filters" snapshot to restore later. Data-independent and portable: genes the
  # other dataset lacks are silently dropped, sliders clamp to their own ranges.

  # Collect the current filter inputs into a base64-encoded JSON code. Empty
  # selections are stored as empty arrays (not dropped) so "show none" round-trips.
  make_filter_code <- function() {
    oe <- function(x) if (is.null(x)) character(0) else x
    vals <- list(
      v             = 1L,
      tab           = input$main_tabs %||% "",
      tier          = oe(input$tier),
      genes         = oe(input$genes),
      impact        = oe(input$impact),
      type          = oe(input$type),
      clnsig        = oe(input$clnsig),
      exclude_am_benign = isTRUE(input$exclude_am_benign),
      cadd          = input$cadd %||% 0,
      revel         = input$revel %||% 0,
      gnomad        = input$gnomad %||% NA,
      min_flags     = input$min_flags %||% 0,
      priority_cadd = input$priority_cadd %||% 20,
      sample_group  = oe(input$sample_group),
      exclude_samples = oe(input$exclude_samples)
    )
    json <- jsonlite::toJSON(vals, auto_unbox = TRUE)
    gsub("[\r\n]", "", jsonlite::base64_enc(charToRaw(as.character(json))))
  }

  # Decode a pasted code and apply each saved filter to this session's inputs.
  apply_filter_code <- function(code) {
    code <- trimws(code %||% "")
    if (!nzchar(code)) {
      showNotification("Paste a filter code first.", type = "warning")
      return(invisible(FALSE))
    }
    vals <- tryCatch(
      jsonlite::fromJSON(rawToChar(jsonlite::base64_dec(code)),
                         simplifyVector = TRUE),
      error = function(e) NULL)
    if (is.null(vals) || is.null(vals$v)) {
      showNotification("That doesn't look like a valid filter code — check it was copied in full.",
                       type = "error", duration = 8)
      return(invisible(FALSE))
    }
    chr <- function(x) as.character(unlist(x))
    if (!is.null(vals$tier))
      updateCheckboxGroupInput(session, "tier", selected = chr(vals$tier))
    if (!is.null(vals$genes))
      updateSelectizeInput(session, "genes", selected = chr(vals$genes))
    if (!is.null(vals$impact))
      updateCheckboxGroupInput(session, "impact", selected = chr(vals$impact))
    if (!is.null(vals$type))
      updateCheckboxGroupInput(session, "type", selected = chr(vals$type))
    if (!is.null(vals$clnsig))
      updateCheckboxGroupInput(session, "clnsig", selected = chr(vals$clnsig))
    updateCheckboxInput(session, "exclude_am_benign",
                        value = isTRUE(vals$exclude_am_benign))
    if (!is.null(vals$cadd))
      updateSliderInput(session, "cadd", value = vals$cadd)
    if (!is.null(vals$revel))
      updateSliderInput(session, "revel", value = vals$revel)
    if (!is.null(vals$gnomad))
      updateNumericInput(session, "gnomad", value = vals$gnomad)
    if (!is.null(vals$min_flags))
      updateSliderInput(session, "min_flags", value = vals$min_flags)
    if (!is.null(vals$priority_cadd))
      updateSliderInput(session, "priority_cadd", value = vals$priority_cadd)
    if (!is.null(vals$sample_group))
      updateCheckboxGroupInput(session, "sample_group", selected = chr(vals$sample_group))
    if (!is.null(vals$exclude_samples))
      updateSelectizeInput(session, "exclude_samples", selected = chr(vals$exclude_samples))
    if (!is.null(vals$tab) && nzchar(vals$tab))
      bslib::nav_select("main_tabs", vals$tab, session = session)
    showNotification("Filters applied from code.", type = "message", duration = 4)
    removeModal()
    invisible(TRUE)
  }

  # Sidebar button: open the share/save dialog with this session's code ready to
  # copy, plus a box to paste someone else's code.
  observeEvent(input$filter_share, {
    code <- make_filter_code()
    showModal(modalDialog(
      title = tagList(bsicons::bs_icon("sliders"), " Share or save filters"),
      tags$p(class = "text-muted",
             "Copy this code to share your current filter settings or save them ",
             "for later. The recipient pastes it below to apply the same filters ",
             "to their own loaded data."),
      tags$label(class = "fw-bold", "Your current filters"),
      tags$div(
        class = "input-group mb-2",
        tags$textarea(id = "filter_code_out", class = "form-control",
                      rows = 2, readonly = NA,
                      style = "font-family:monospace;font-size:0.8rem;",
                      code),
        tags$button(
          tagList(bsicons::bs_icon("clipboard"), " Copy"),
          class = "btn btn-outline-primary",
          onclick = paste0(
            "var t=document.getElementById('filter_code_out');",
            "t.select();navigator.clipboard.writeText(t.value);",
            "this.innerHTML='Copied!';return false;"))
      ),
      tags$hr(),
      textAreaInput("filter_code_in", tags$span(class = "fw-bold", "Apply a filter code"),
                    value = "", rows = 2, width = "100%",
                    placeholder = "Paste a filter code here…"),
      easyClose = TRUE,
      footer = tagList(
        actionButton("filter_code_apply",
                     tagList(bsicons::bs_icon("check2-circle"), " Apply filters"),
                     class = "btn btn-primary"),
        modalButton("Close"))
    ))
  })

  # Apply button inside the dialog decodes and applies the pasted code.
  observeEvent(input$filter_code_apply, {
    apply_filter_code(input$filter_code_in)
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

  # Action / external-resource buttons shown with the active variant:
  #  * "View all variants" clears the selected variant and returns to the
  #    whole-gene protein lollipop landing (with its variant table).
  #  * "ClinVar record" opens the variation page in the real default browser
  #    (server-side handler, so it doesn't get trapped in RStudio's Viewer).
  # Re-renders with the active variant so the ClinVar link tracks the selection.
  output$variant_links <- renderUI({
    row <- modal_variant(); req(row)
    open_btn <- function(url, label, icon) tags$button(
      tagList(bsicons::bs_icon(icon), paste0(" ", label)),
      type = "button", class = "btn btn-sm btn-outline-primary mb-2 me-2",
      onclick = sprintf(
        "Shiny.setInputValue('open_url','%s',{priority:'event'});", .jsesc(url)))

    # Return to the whole-gene lollipop landing (deselect the current variant).
    view_all_btn <- tags$button(
      tagList(bsicons::bs_icon("graph-up"), " View all variants"),
      type = "button", class = "btn btn-sm btn-outline-secondary mb-2 me-2",
      onclick = "Shiny.setInputValue('view_all_variants',Math.random(),{priority:'event'});")

    # ClinVar: deep-link to the variation record when a CLNVID is present.
    cv_id <- if ("CLNVID" %in% names(row))
      suppressWarnings(as.character(row$CLNVID)) else NA
    cv_id <- if (length(cv_id) && !is.na(cv_id) && nzchar(cv_id) &&
                 cv_id != "NA") cv_id else NA
    cv_btn <- if (!is.na(cv_id)) open_btn(
      sprintf("https://www.ncbi.nlm.nih.gov/clinvar/variation/%s/", cv_id),
      sprintf("ClinVar record (%s)", cv_id), "clipboard2-pulse") else NULL

    div(class = "mb-1", view_all_btn, cv_btn)
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

    chips <- lapply(carriers, function(fid) {
      tags$a(href = "#", fmt_sample(fid),
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
          legend_dot("#C62828", "MacTel"), legend_dot("#1565C0", "HSAN1"),
          legend_dot("#6A1B9A", "MacTel + HSAN1"),   legend_dot("#546E7A", "Control")),
      div(class = "mb-2", chips)
    )
  })

  # ---- 3D AlphaFold structure viewer ----------------------------------------
  # UI block dropped into the gene/variant modals. Returns NULL when r3dmol is
  # not installed, so the rest of the modal still renders.
  structure_section <- function() {
    if (!HAS_R3DMOL) return(NULL)
    tagList(
      tags$hr(class = "mt-3"),
      tags$div(
        class = "mb-3",
        tags$strong("AlphaFold structure"),
        tags$span(class = "text-muted small",
                  " — cartoon coloured by per-residue confidence (pLDDT); ",
                  "the selected variant residue is highlighted in magenta.")),
      r3dmol::r3dmolOutput("variant_structure", height = "440px")
    )
  }

  # Renders the AlphaFold model for the active gene, coloured by pLDDT, with the
  # selected variant residue drawn as magenta sticks + sphere and zoomed-to.
  # Re-renders whenever the active gene or variant changes (e.g. clicking a
  # different lollipop point). Reused by both the gene and variant modals.
  if (HAS_R3DMOL) output$variant_structure <- r3dmol::renderR3dmol({
    gene <- modal_gene(); req(gene)
    uni  <- af_uniprot_for_gene(gene, PROTEIN_DOMAINS)
    validate(need(!is.na(uni),
                  "No UniProt accession is mapped for this gene, so no structure is available."))
    pdb <- fetch_alphafold_pdb(uni)
    validate(need(!is.null(pdb),
                  "AlphaFold structure could not be loaded (offline, or this protein is not modelled)."))

    row  <- modal_variant()
    resi <- if (!is.null(row) && nrow(row) > 0)
      aa_position(row$HGVSp_short) else NA_integer_

    v <- r3dmol::r3dmol(viewer_spec = r3dmol::m_viewer_spec(
           cartoonQuality = 8, backgroundColor = "#FFFFFF")) %>%
      r3dmol::m_add_model(data = pdb, format = "pdb") %>%
      r3dmol::m_set_style(style = r3dmol::m_style_cartoon(
        colorfunc = htmlwidgets::JS(AF_PLDDT_COLORFUNC)))

    if (!is.na(resi)) {
      # "magentaCarbon" colours carbons magenta (heteroatoms keep element
      # colours) so the variant residue pops against the cartoon. r3dmol always
      # emits colorscheme:"default", which 3Dmol prioritises over a plain
      # `color`, so we set the scheme rather than a single colour. Spheres are
      # scaled up so the residue is easy to spot in the whole-protein view.
      sel <- r3dmol::m_sel(resi = resi)
      v <- v %>%
        r3dmol::m_add_style(sel = sel,
          style = r3dmol::m_style_stick(colorScheme = "magentaCarbon", radius = 0.4)) %>%
        r3dmol::m_add_style(sel = sel,
          style = r3dmol::m_style_sphere(colorScheme = "magentaCarbon", scale = 0.9))
    }
    # Always frame the entire protein (never zoom to the residue) so the view is
    # identical for every variant in the gene — jumping between variants then
    # only moves the magenta marker, making the position change easy to follow.
    # A fixed zoom factor on the whole-protein framing pulls the camera in a
    # little (so the structure fills more of the viewer) while staying identical
    # across variants.
    v %>%
      r3dmol::m_zoom_to(sel = r3dmol::m_sel()) %>%
      r3dmol::m_zoom(factor = 1.4)
  })

  # One protein modal with exactly two states, driven by modal_variant():
  #   * gene state    (row = NULL): clickable variant table above the lollipop,
  #                                 no ring, straight gene-report download.
  #   * variant state (row set):    variant detail + carriers + "View all
  #                                 variants" button, ringed lollipop, report
  #                                 chooser.
  # Every variant click (table row OR lollipop point) just sets modal_variant(),
  # so the page toggles between these two states in place — there is no third
  # hybrid state, and clicking via the plot or the table does the same thing.
  show_protein_modal <- function(gene, row = NULL) {
    if (is.null(gene) || is.na(gene) || gene == "") return(invisible())
    modal_gene(gene)
    modal_variant(if (!is.null(row) && nrow(row) > 0) row else NULL)

    gi <- gene_info_rv()
    info <- if (!is.null(gi)) gi[gi$SYMBOL == gene, ] else NULL
    tier <- if (!is.null(info) && nrow(info) > 0) info$Tier[1]
            else if (!is.null(row) && !is.null(row$Tier)) row$Tier else NA
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
        tags$span(gene, style = "font-weight:700;font-size:1.2rem;"),
        if (!is.null(tier) && !is.na(tier))
          tags$span(tier, class = "badge bg-secondary",
                    style = "margin-left:8px;vertical-align:middle;")
      ),
      uiOutput("variant_detail"),
      gene_desc,
      uiOutput("variant_links"),
      uiOutput("variant_samples"),
      # Gene state only: the clickable variant table above the lollipop.
      conditionalPanel(
        condition = "output.modal_gene_state === true",
        tags$p(tags$strong("Variants in this gene"),
               style = "margin-bottom:4px;"),
        uiOutput("legend_gene_variants"),
        DT::DTOutput("gene_variant_table")
      ),
      plotly::plotlyOutput("lollipop", height = 430),
      div(class = "d-flex justify-content-end mt-1",
          downloadButton("dl_lollipop", "Download plot (PNG)",
                         class = "btn-sm btn-outline-secondary",
                         icon = bsicons::bs_icon("image"))),
      structure_section(),
      easyClose = TRUE, size = "xl",
      footer = uiOutput("report_footer")
    ))
  }

  # JS-visible flag for the conditionalPanel: TRUE in gene state (no variant).
  output$modal_gene_state <- reactive({ is.null(modal_variant()) })
  outputOptions(output, "modal_gene_state", suspendWhenHidden = FALSE)

  # Footer swaps with the state: gene state downloads the gene report directly;
  # variant state asks which report (variant vs gene) via report_choose.
  output$report_footer <- renderUI({
    dl_icon <- bsicons::bs_icon("file-earmark-arrow-down")
    if (is.null(modal_variant())) {
      tagList(
        downloadButton("report_dl", tagList(dl_icon, " Download report"),
                       class = "btn btn-primary", icon = NULL),
        modalButton("Close"))
    } else {
      tagList(
        actionButton("report_choose", tagList(dl_icon, " Download report"),
                     class = "btn btn-primary"),
        modalButton("Close"))
    }
  })

  output$legend_gene_variants <- renderUI({
    gene <- modal_gene(); req(gene)
    diag_legend(dplyr::filter(raw(), SYMBOL == gene)$family_id)
  })

  # Clickable variant table shown on the gene-landing modal (above the lollipop).
  # Six columns matching the gene report; Variant and Sample cells are clickable
  # (Gene is omitted since the whole modal is already scoped to one gene).
  output$gene_variant_table <- DT::renderDT({
    gene <- modal_gene(); req(gene)
    gdf <- dplyr::filter(raw(), SYMBOL == gene)
    validate(need(nrow(gdf) > 0, "No variants for this gene."))
    tbl <- gdf %>%
      dplyr::group_by(CHROM, POS, REF, ALT) %>%
      dplyr::summarise(
        HGVSp    = dplyr::first(HGVSp_short),
        Impact   = as.character(dplyr::first(IMPACT)),
        CADD     = suppressWarnings(max(CADD, na.rm = TRUE)),
        ClinVar  = as.character(dplyr::first(CLNSIG_clean)),
        carriers = list(sort(unique(family_id))),
        .groups  = "drop") %>%
      dplyr::arrange(dplyr::desc(CADD))
    sample_links <- vapply(tbl$carriers,
      function(fs) paste(link_sample(fs), collapse = " "), character(1))
    out <- data.frame(
      Variant = link_variant(tbl$CHROM, tbl$POS, tbl$REF, tbl$ALT,
                             input_id = "modal_pick_variant"),
      HGVSp   = tbl$HGVSp,
      Impact  = tbl$Impact,
      CADD    = ifelse(is.finite(tbl$CADD), round(tbl$CADD, 1), NA_real_),
      ClinVar = tbl$ClinVar,
      Samples = sample_links,
      check.names = FALSE, stringsAsFactors = FALSE)
    names(out)[names(out) == "Samples"] <- "Sample(s)"
    DT::datatable(out, escape = FALSE, rownames = FALSE, selection = "none",
                  filter = "top",
                  options = list(pageLength = 10, scrollX = TRUE,
                                 order = order_desc_by(out, "CADD"),
                                 dom = "ftip")) %>%
      DT::formatStyle("Impact",
                      backgroundColor = DT::styleEqual("HIGH", "#FDDEDE")) %>%
      DT::formatStyle("ClinVar",
                      backgroundColor = DT::styleEqual(
                        c("Pathogenic", "Pathogenic/Likely_pathogenic",
                          "Likely_pathogenic"),
                        c("#FDDEDE", "#F7D6F7", "#FDE8D8")))
  })

  # ---- one-click report (self-contained HTML) ------------------------------
  # Build a standalone, emailable HTML summary for the gene/variant currently
  # shown in the protein modal. The lollipop is embedded as a base64 PNG so the
  # file is fully self-contained — no internet or extra files needed to view it.
  make_report_html <- function(gene, row) {
    esc <- function(x) htmltools::htmlEscape(as.character(x))
    has_row <- !is.null(row) && nrow(row) > 0

    gi <- gene_info_rv()
    ginfo <- if (!is.null(gi)) gi[gi$SYMBOL == gene, ] else NULL
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
    # Gene report (no single variant): label every position on the lollipop.
    p <- tryCatch(plot_variant_lollipop(gdf, ddf, gene, sel_key,
                                        label_all = !has_row,
                                        italic_gene = TRUE),
                  error = function(e) NULL)
    plot_html <- "<p class='muted'>No protein-coding positions to plot for this gene.</p>"
    if (!is.null(p)) {
      # For the static report, put the legend below the plot so the lollipop
      # uses the full page width (the right-hand legend squeezes it otherwise).
      p_png <- p +
        ggplot2::theme(
          legend.position = "bottom",
          legend.box      = "vertical",
          legend.title    = ggplot2::element_text(size = 9),
          legend.text     = ggplot2::element_text(size = 8)) +
        ggplot2::guides(
          colour = ggplot2::guide_legend(order = 1, nrow = 2, byrow = TRUE,
                                         override.aes = list(size = 3.5)),
          size   = ggplot2::guide_legend(order = 2, nrow = 1),
          fill   = ggplot2::guide_legend(order = 3, ncol = 1))
      tmp <- tempfile(fileext = ".png")
      ok <- tryCatch({
        ggplot2::ggsave(tmp, p_png, width = 9, height = 5.6, dpi = 130, bg = "white")
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
        if (is.null(SAMPLE_INFO)) return("#546E7A")
        r <- SAMPLE_INFO[SAMPLE_INFO$family_id == fid, , drop = FALSE]
        if (nrow(r) == 0) return("#546E7A")
        m <- isTRUE(r$is_mactel[1]); h <- isTRUE(r$is_hsan1[1])
        if (m && h) "#6A1B9A" else if (m) "#C62828" else if (h) "#1565C0" else "#546E7A"
      }
      if (length(carriers) > 0) {
        pills <- paste(vapply(carriers, function(fid)
          sprintf('<span class="pill" style="background:%s">%s</span>',
                  diag_colour(fid), esc(fmt_sample(fid))), character(1)), collapse = "")
        body_html <- paste0(body_html,
          sprintf("<h3>Carried by %d sample%s</h3>", length(carriers),
                  if (length(carriers) == 1) "" else "s"),
          "<p class='legend'>",
          "<span class='dot' style='background:#C62828'></span>MacTel",
          "<span class='dot' style='background:#1565C0'></span>HSAN1",
          "<span class='dot' style='background:#6A1B9A'></span>MacTel + HSAN1",
          "<span class='dot' style='background:#546E7A'></span>Control</p>",
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
                         # one row per variant; list every sample carrying it
                         Samples = paste(sort(unique(fmt_sample(family_id))), collapse = ", "),
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
        "<h3>%d variant%s in this gene</h3><table class='vt'><thead><tr><th>Variant</th><th>HGVSp</th><th>Impact</th><th>CADD</th><th>ClinVar</th><th>Sample(s)</th></tr></thead><tbody>%s</tbody></table>",
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

  # Gene-report filename (whole gene, no single variant).
  gene_report_name <- function() {
    g <- modal_gene(); g <- if (is.null(g)) "gene" else g
    sprintf("MacTel_report_%s_%s.html", g, Sys.Date())
  }
  # Variant-report filename (one selected variant).
  variant_report_name <- function() {
    g <- modal_gene(); r <- modal_variant()
    g <- if (is.null(g)) "gene" else g
    if (!is.null(r) && nrow(r) > 0)
      sprintf("MacTel_report_%s_%s-%s_%s.html", g, r$CHROM, r$POS, Sys.Date())
    else gene_report_name()
  }
  write_report <- function(file, row) {
    g <- modal_gene()
    if (is.null(g)) {
      writeLines("<html><body><p>No gene selected.</p></body></html>", file)
      return(invisible())
    }
    writeLines(make_report_html(g, row), file)
  }

  # Gene-landing modal footer: straight gene-report download (no variant chosen).
  output$report_dl <- downloadHandler(
    filename = gene_report_name,
    content  = function(file) write_report(file, modal_variant())
  )

  # Variant modal "Download report" -> ask whether they want the single-variant
  # report or the whole-gene report (with every position labelled on the plot).
  observeEvent(input$report_choose, {
    req(modal_variant())
    showModal(modalDialog(
      title = tagList(bsicons::bs_icon("file-earmark-arrow-down"),
                      " Download report"),
      tags$p("Which report would you like?"),
      tags$ul(
        tags$li(tags$strong("Variant report"),
                " — details, carriers and protein lollipop for the selected variant."),
        tags$li(tags$strong("Gene report"),
                " — every variant and its samples for this gene, plus the full ",
                "protein lollipop with each position labelled.")),
      easyClose = TRUE, size = "m",
      footer = tagList(
        downloadButton("report_dl_variant",
                       tagList(bsicons::bs_icon("crosshair"), " Variant report"),
                       class = "btn btn-primary"),
        downloadButton("report_dl_gene",
                       tagList(bsicons::bs_icon("diagram-3"), " Gene report"),
                       class = "btn btn-outline-primary"),
        modalButton("Cancel"))
    ))
  })

  output$report_dl_variant <- downloadHandler(
    filename = variant_report_name,
    content  = function(file) { write_report(file, modal_variant()); removeModal() }
  )
  output$report_dl_gene <- downloadHandler(
    filename = gene_report_name,
    content  = function(file) { write_report(file, NULL); removeModal() }
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
    gg <- plotly::ggplotly(p, tooltip = "text", source = "lollipop") %>%
      plotly::event_register("plotly_click") %>%
      # ggplotly cannot render a plotmath italic title, so italicise the gene
      # symbol here with an HTML tag (the subtitle is dropped by ggplotly anyway).
      plotly::layout(title = list(
        text = sprintf("<i>%s</i> protein lollipop", gene)))
    # ggplotly drops ggplot annotations, so re-add the non-coding disclaimer
    # as a native plotly annotation when the selected variant can't be drawn.
    if (isTRUE(attr(p, "sel_not_coding"))) {
      gg <- gg %>% plotly::layout(annotations = list(list(
        x = 0.5, y = 1, xref = "paper", yref = "paper",
        xanchor = "center", yanchor = "top",
        text = "Selected variant is not protein coding",
        showarrow = FALSE,
        font = list(color = "#664d03", size = 14),
        bgcolor = "#fff3cd", bordercolor = "#664d03",
        borderwidth = 1, borderpad = 4)))
    }
    gg
  })

  # Static PNG of the protein lollipop for the gene currently in the modal.
  # Every position is labelled (as in the gene report) so the image stands
  # alone without the interactive hover, and the legend sits below the plot so
  # the lollipop uses the full width.
  output$dl_lollipop <- downloadHandler(
    filename = function() sprintf("%s_lollipop_%s.png",
                                  modal_gene() %||% "gene", Sys.Date()),
    content  = function(file) {
      gene <- modal_gene(); req(gene)
      gdf  <- dplyr::filter(raw(), SYMBOL == gene)
      ddf  <- if (!is.null(PROTEIN_DOMAINS))
        dplyr::filter(PROTEIN_DOMAINS, SYMBOL == gene) else NULL
      row     <- modal_variant()
      sel_key <- if (!is.null(row) && nrow(row) > 0)
        paste(row$CHROM, row$POS, row$REF, row$ALT) else NULL
      p <- tryCatch(plot_variant_lollipop(gdf, ddf, gene, sel_key,
                                          label_all = TRUE,
                                          italic_gene = TRUE),
                    error = function(e) NULL)
      req(!is.null(p))
      p <- p +
        ggplot2::theme(
          legend.position = "bottom",
          legend.box      = "vertical",
          legend.title    = ggplot2::element_text(size = 9),
          legend.text     = ggplot2::element_text(size = 8)) +
        ggplot2::guides(
          colour = ggplot2::guide_legend(order = 1, nrow = 2, byrow = TRUE,
                                         override.aes = list(size = 3.5)),
          size   = ggplot2::guide_legend(order = 2, nrow = 1),
          fill   = ggplot2::guide_legend(order = 3, ncol = 1))
      ggplot2::ggsave(file, p, device = "png", width = 11, height = 6.4,
                      dpi = 200, bg = "white")
    }
  )

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
    show_gene_modal(input$cell_gene, gene_info_rv())
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

  # "View variants on protein" in the gene modal -> open the protein modal in
  # GENE state (variant table above the lollipop, no variant selected).
  observeEvent(input$gene_view_lollipop, {
    show_protein_modal(input$gene_view_lollipop, NULL)
  })

  # "View all variants" in the variant state -> switch the already-open modal
  # back to GENE state in place (no re-open, smooth toggle).
  observeEvent(input$view_all_variants, {
    req(modal_gene())
    modal_variant(NULL)
  })

  # Clicking a Variant cell in the whole-gene table (inside the modal) -> switch
  # the already-open modal to VARIANT state in place (no re-open, smooth toggle).
  observeEvent(input$modal_pick_variant, {
    g <- modal_gene(); req(g)
    parts <- strsplit(input$modal_pick_variant, "\\|\\|")[[1]]
    if (length(parts) != 4) return()
    hit <- raw() %>%
      dplyr::filter(SYMBOL == g, CHROM == parts[1], as.character(POS) == parts[2],
                    REF == parts[3], ALT == parts[4])
    if (nrow(hit) > 0) modal_variant(hit[1, ])
  })

  # Clicking a Variant cell anywhere OUTSIDE the modal (main variant table,
  # priority table, etc.) -> open the protein modal in VARIANT state.
  observeEvent(input$cell_variant, {
    parts <- strsplit(input$cell_variant, "\\|\\|")[[1]]
    if (length(parts) != 4) return()
    hit <- raw() %>%
      dplyr::filter(CHROM == parts[1], as.character(POS) == parts[2],
                    REF == parts[3], ALT == parts[4])
    if (nrow(hit) > 0) show_protein_modal(hit$SYMBOL[1], hit[1, ])
  })

  # Clicking a Sample cell (table or variant-modal chip) -> close any modal and
  # jump to the Sample explorer for that sample.
  observeEvent(input$cell_sample, {
    removeModal()
    fid <- input$cell_sample
    if (is.null(fid) || !nzchar(fid)) return()
    # The picker only lists samples that carry candidate variants. A family
    # member with none (still clickable from the Family explorer header) is not
    # among those choices, so updateSelectizeInput(selected = fid) would
    # silently clear the selection and land on an empty Sample explorer. Add the
    # fid to the choices first so the sample (with its identity + empty variant
    # table) loads correctly.
    df <- raw()
    fids <- if (is.null(df)) fid else sort(unique(c(df$family_id, fid)))
    updateSelectizeInput(session, "sample_pick",
                         choices = stats::setNames(fids, fmt_sample(fids)),
                         selected = fid, server = TRUE)
    bslib::nav_select("main_tabs", "Sample explorer")
  })

  # Clicking a Family_ID badge -> remember the family, sync the dropdown, and
  # jump to the Family explorer tab.
  observeEvent(input$cell_family, {
    removeModal()
    fam <- input$cell_family
    if (is.null(fam) || !nzchar(fam)) return()
    selected_family(fam)
    updateSelectizeInput(session, "family_pick", selected = fam)
    bslib::nav_select("main_tabs", "Family explorer")
  })

  # Picking a family straight from the Family explorer dropdown (no click needed
  # from the Sample explorer). Empty selection clears back to the placeholder.
  observeEvent(input$family_pick, {
    fam <- input$family_pick
    selected_family(if (is.null(fam) || !nzchar(fam)) NULL else fam)
  }, ignoreNULL = FALSE, ignoreInit = TRUE)

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
        dplyr::mutate(family_id = fmt_sample(family_id)) %>%
        dplyr::rename(Sample = family_id)
      readr::write_csv(d, file)
    }
  )

  # ---- Excel export with column picker --------------------------------------
  # Both the Variant table and Priority variants tabs share one modal: a column
  # picker (offering every column in the loaded data, with a sensible default
  # ticked) that produces a nicely formatted .xlsx.
  export_target <- reactiveVal("variants")     # "variants" | "priority"

  # Lead-column order mirroring the interactive variant table (display name ->
  # underlying column). Any columns not listed here keep their original order
  # after these.
  EXPORT_LEAD_ORDER <- c("SYMBOL", "Tier", "family_id", "Variant", "HGVSc",
                         "HGVSp", "HGVSp_short", "Consequence", "IMPACT",
                         "TYPE", "CADD", "REVEL", "am_class", "SpliceAI_max",
                         "CLNSIG_clean", "gnomad_AF", "inheritance")

  export_dataset <- function(target) {
    d <- if (identical(target, "priority")) priority() else filtered()
    if (is.null(d) || nrow(d) == 0) return(d)
    # Synthesise the same "Variant" column the interactive table shows
    # (CHROM:POS REF>ALT), order the leading columns to match that table, and
    # rank rows by descending CADD (NAs last) as everywhere else in the app.
    d %>%
      dplyr::mutate(Variant = sprintf("%s:%s %s>%s", CHROM, POS, REF, ALT)) %>%
      dplyr::relocate(dplyr::any_of(EXPORT_LEAD_ORDER)) %>%
      dplyr::arrange(dplyr::desc(CADD))
  }

  export_default <- function(target, cols) {
    def <- EXPORT_DEFAULT_COLS
    if (identical(target, "priority")) def <- c(def, EXPORT_PRIORITY_EXTRA)
    intersect(def, cols)
  }

  show_export_modal <- function(target) {
    d <- export_dataset(target)
    if (is.null(d) || nrow(d) == 0) {
      showModal(modalDialog(
        title = "Nothing to export", easyClose = TRUE,
        "No variants match the current filters.",
        footer = modalButton("Close")))
      return(invisible())
    }
    cols    <- names(d)
    labels  <- cols
    labels[labels == "family_id"]   <- "Sample (family_id)"
    labels[labels == "HGVSp_short"] <- "HGVSp (protein change, no prefix)"
    choices <- stats::setNames(cols, labels)
    label   <- if (identical(target, "priority")) "priority variants" else
               "filtered variants"
    showModal(modalDialog(
      title = sprintf("Export %s to Excel", label), size = "l",
      easyClose = FALSE,
      tags$p(class = "text-muted small",
             sprintf("%s rows. Tick the columns to include, then download a formatted .xlsx table.",
                     format(nrow(d), big.mark = ","))),
      tags$style(HTML(".export-cols .checkbox{break-inside:avoid;margin:0 0 2px;}")),
      div(class = "mb-2",
          actionButton("export_all", "Select all",
                       class = "btn-sm btn-outline-secondary"),
          actionButton("export_none", "Clear all",
                       class = "btn-sm btn-outline-secondary ms-1"),
          actionButton("export_reset", "Reset to default",
                       class = "btn-sm btn-outline-secondary ms-1")),
      div(class = "export-cols",
          style = paste("column-count:3;column-gap:1.5rem;max-height:55vh;",
                        "overflow-y:auto;border:1px solid #dee2e6;",
                        "border-radius:6px;padding:10px;"),
          checkboxGroupInput("export_cols", NULL, choices = choices,
                             selected = export_default(target, cols),
                             width = "100%")),
      footer = tagList(
        modalButton("Cancel"),
        downloadButton("do_export_xlsx", "Download .xlsx",
                       class = "btn-success"))
    ))
  }

  observeEvent(input$xl_table,
               { export_target("variants"); show_export_modal("variants") })
  observeEvent(input$xl_priority,
               { export_target("priority"); show_export_modal("priority") })

  observeEvent(input$export_all, {
    d <- export_dataset(export_target()); req(d)
    updateCheckboxGroupInput(session, "export_cols", selected = names(d))
  })
  observeEvent(input$export_none, {
    updateCheckboxGroupInput(session, "export_cols", selected = character(0))
  })
  observeEvent(input$export_reset, {
    d <- export_dataset(export_target()); req(d)
    updateCheckboxGroupInput(session, "export_cols",
                             selected = export_default(export_target(), names(d)))
  })

  output$do_export_xlsx <- downloadHandler(
    filename = function() sprintf("%s_%s.xlsx",
      if (identical(export_target(), "priority")) "priority_variants" else
        "filtered_variants", Sys.Date()),
    content = function(file) {
      target <- export_target()
      d <- export_dataset(target); req(d)
      sel <- input$export_cols
      if (is.null(sel) || length(sel) == 0)
        sel <- export_default(target, names(d))
      sel <- intersect(names(d), sel)          # keep original column order
      out <- d[, sel, drop = FALSE]
      if ("family_id" %in% names(out)) {
        out$family_id <- fmt_sample(out$family_id)
        names(out)[names(out) == "family_id"] <- "Sample"
      }
      # HGVSp_short is the protein change without the transcript prefix (as in
      # the gene report). Give it the clean "HGVSp" header, unless the raw
      # prefixed HGVSp is also included (then keep the names distinct).
      if ("HGVSp_short" %in% names(out)) {
        new_name <- if ("HGVSp" %in% names(out)) "HGVSp (no prefix)" else "HGVSp"
        names(out)[names(out) == "HGVSp_short"] <- new_name
      }
      out[] <- lapply(out, function(x) if (is.factor(x)) as.character(x) else x)
      write_variants_xlsx(out, file, sheet =
        if (identical(target, "priority")) "Priority variants" else "Variants")
      removeModal()
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
        dplyr::mutate(family_id = fmt_sample(family_id)) %>%
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

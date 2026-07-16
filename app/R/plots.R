# =============================================================================
# R/plots.R — reusable ggplot builders for the variant explorer
# Each function takes an already-filtered dataframe and returns a ggplot.
# =============================================================================

theme_app <- function(...) {
  ggplot2::theme_bw(base_size = 12) +
    ggplot2::theme(
      plot.title      = ggplot2::element_text(face = "bold", size = 13),
      plot.subtitle   = ggplot2::element_text(size = 10, colour = "grey40"),
      axis.title      = ggplot2::element_text(face = "bold"),
      legend.position = "right"
    ) +
    ggplot2::theme(...)
}

# Angled x-axis tick labels, shared by the categorical overview bar charts so
# crowded category names do not overlap.
angle_x <- function() ggplot2::element_text(angle = 45, hjust = 1)

# Okabe-Ito qualitative colour-blind-safe palette (8 hues).
OKABE_ITO <- c("#E69F00", "#56B4E9", "#009E73", "#F0E442",
               "#0072B2", "#D55E00", "#CC79A7", "#000000")

# Remap a named default palette to an alternative colour scheme while keeping
# the same category names/levels. "Default" returns the semantic palette
# unchanged; "Colour-blind" uses the Okabe-Ito set; any other value is treated
# as a base-R grDevices::hcl.colors palette name (e.g. "Viridis", "Cividis",
# "Set 2"). Unknown names fall back to the semantic palette. No extra package
# dependencies.
apply_palette <- function(default_named, palette = "Default") {
  if (is.null(palette) || palette == "Default") return(default_named)
  n <- length(default_named)
  cols <- if (palette == "Colour-blind") {
    rep(OKABE_ITO, length.out = n)
  } else {
    tryCatch(grDevices::hcl.colors(n, palette),
             error = function(e) unname(default_named))
  }
  stats::setNames(unname(cols), names(default_named))
}

# Return a count-axis scale for a categorical count bar chart. When the counts
# span a wide range (largest category >= 20x the smallest non-zero one), a
# linear axis buries the small-but-informative bars, so switch to a base-10
# pseudo-log scale (pseudo-log keeps counts of 0 and 1 sensible, unlike plain
# log10). Otherwise use a linear axis. `axis` selects the count axis: "y" for
# vertical bars, "x" for horizontal bars.
dynamic_count_y_scale <- function(n,
                                  expand = ggplot2::expansion(mult = c(0, 0.18)),
                                  axis = c("y", "x")) {
  axis     <- match.arg(axis)
  nz       <- n[n > 0]
  use_log  <- length(nz) >= 2 && (max(nz) / min(nz)) >= 20
  scale_fn <- if (axis == "x") ggplot2::scale_x_continuous
              else             ggplot2::scale_y_continuous
  if (use_log) {
    scale_fn(
      transform = scales::pseudo_log_trans(base = 10),
      breaks    = c(0, 1, 3, 10, 30, 100, 300, 1000),
      expand    = expand)
  } else {
    scale_fn(expand = expand)
  }
}

# --- IMPACT distribution -----------------------------------------------------
plot_impact <- function(df, palette = "Default") {
  cnt <- dplyr::count(df, IMPACT, .drop = FALSE)
  ggplot2::ggplot(cnt, ggplot2::aes(IMPACT, n, fill = IMPACT)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.4,
                       fontface = "bold", size = 3.5) +
    ggplot2::scale_fill_manual(values = apply_palette(COL_IMPACT, palette),
                               drop = FALSE) +
    dynamic_count_y_scale(cnt$n) +
    ggplot2::labs(title = "VEP impact", x = NULL, y = "Variants") +
    theme_app(legend.position = "none", axis.text.x = angle_x())
}

# --- Variant TYPE ------------------------------------------------------------
plot_type <- function(df, palette = "Default") {
  cnt <- dplyr::count(df, TYPE, .drop = FALSE)
  ggplot2::ggplot(cnt, ggplot2::aes(TYPE, n, fill = TYPE)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.4,
                       fontface = "bold", size = 3.5) +
    ggplot2::scale_fill_manual(values = apply_palette(COL_TYPE, palette),
                               drop = FALSE) +
    dynamic_count_y_scale(cnt$n) +
    ggplot2::labs(title = "Variant type", x = NULL, y = "Variants") +
    theme_app(legend.position = "none", axis.text.x = angle_x())
}

# --- ClinVar -----------------------------------------------------------------
plot_clnsig <- function(df, palette = "Default") {
  cnt <- dplyr::count(df, CLNSIG_clean, .drop = FALSE)

  ggplot2::ggplot(cnt, ggplot2::aes(CLNSIG_clean, n, fill = CLNSIG_clean)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.4,
                       fontface = "bold", size = 3.2) +
    ggplot2::scale_fill_manual(values = apply_palette(COL_CLNSIG, palette),
                               drop = FALSE) +
    dynamic_count_y_scale(cnt$n) +
    ggplot2::scale_x_discrete(labels = function(x)
      stringr::str_wrap(gsub("_", " ", x), 12)) +
    ggplot2::labs(title = "ClinVar classification", x = NULL, y = "Variants") +
    theme_app(legend.position = "none",
              axis.text.x = ggplot2::element_text(angle = 45, hjust = 1,
                                                  size = 8))
}

# --- CADD histogram ----------------------------------------------------------
plot_cadd <- function(df, threshold = 20, palette = "Default") {
  d <- dplyr::filter(df, !is.na(CADD))
  if (nrow(d) == 0) return(NULL)

  # Drop the legend into whichever top corner has the most headroom: bin the
  # CADD values the same way the histogram does (binwidth 2), then compare the
  # tallest bar in the left vs right half of the range. Ties favour the right.
  rng  <- range(d$CADD)
  brks <- seq(floor(rng[1] / 2) * 2, ceiling(rng[2] / 2) * 2, by = 2)
  if (length(brks) > 1) {
    cnts      <- graphics::hist(d$CADD, breaks = brks,
                                include.lowest = TRUE, plot = FALSE)$counts
    mids      <- utils::head(brks, -1) + 1
    mid_x     <- mean(rng)
    left_max  <- max(c(0, cnts[mids <  mid_x]))
    right_max <- max(c(0, cnts[mids >= mid_x]))
  } else {
    left_max <- right_max <- 0
  }
  on_right <- right_max <= left_max
  leg_pos  <- if (on_right) c(0.98, 0.98) else c(0.02, 0.98)
  leg_just <- if (on_right) c(1, 1)       else c(0, 1)

  ggplot2::ggplot(d, ggplot2::aes(CADD, fill = IMPACT)) +
    ggplot2::geom_histogram(binwidth = 2, colour = "white", alpha = 0.9) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed",
                        colour = "red", linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = apply_palette(COL_IMPACT, palette),
                               drop = TRUE, name = "Impact") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(title = "CADD distribution", x = NULL, y = "Count") +
    theme_app(
      legend.position        = "inside",
      legend.position.inside = leg_pos,   # adaptive: emptier top corner
      legend.justification   = leg_just,
      legend.background      = ggplot2::element_rect(
        fill = scales::alpha("white", 0.6), colour = NA),
      legend.key.size        = ggplot2::unit(11, "pt"),   # compact for this panel
      legend.title           = ggplot2::element_text(size = 9),
      legend.text            = ggplot2::element_text(size = 8))
}

# --- Inheritance -------------------------------------------------------------
plot_inheritance <- function(df, palette = "Default") {
  cnt <- dplyr::count(df, inheritance)
  p <- ggplot2::ggplot(cnt, ggplot2::aes(stats::reorder(inheritance, n), n,
                                         fill = inheritance)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n), hjust = -0.2,
                       fontface = "bold", size = 3.5) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Inheritance mode", x = NULL, y = "Variants") +
    theme_app(legend.position = "none")
  # Only override the default hue palette when an alternative is requested.
  if (!is.null(palette) && palette != "Default") {
    lv <- sort(unique(as.character(cnt$inheritance)))
    def <- stats::setNames(seq_along(lv), lv)   # names carry the levels
    p <- p + ggplot2::scale_fill_manual(values = apply_palette(def, palette))
  }
  p
}

# --- VEP consequence ---------------------------------------------------------
#' Horizontal bar chart of the most common VEP consequences. VEP reports one or
#' more &-joined terms per variant (most severe first); we keep the lead term so
#' compound calls collapse to a single, readable category. The top n_top
#' consequences are shown, with any remainder pooled into "other".
plot_consequence <- function(df, n_top = 12, palette = "Default") {
  if (!("Consequence" %in% names(df))) return(NULL)
  d <- df %>%
    dplyr::filter(!is.na(Consequence), Consequence != "") %>%
    dplyr::mutate(cons = sub("&.*", "", as.character(Consequence)))
  if (nrow(d) == 0) return(NULL)

  cnt <- d %>%
    dplyr::count(cons, name = "n") %>%
    dplyr::arrange(dplyr::desc(n))
  if (nrow(cnt) > n_top) {
    cnt <- rbind(cnt[seq_len(n_top), ],
                 data.frame(cons = "other",
                            n = sum(cnt$n[(n_top + 1):nrow(cnt)])))
  }
  cnt$label <- gsub("_", " ", cnt$cons)
  cnt$label <- sub(" variant$", "", cnt$label)          # trim redundant suffix
  cnt$label <- stringr::str_wrap(cnt$label, 20)         # keep the column narrow
  cnt$label <- factor(cnt$label, levels = cnt$label[order(cnt$n)])  # asc y axis

  p <- ggplot2::ggplot(cnt, ggplot2::aes(n, label, fill = label)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n), hjust = -0.2,
                       fontface = "bold", size = 3.5) +
    dynamic_count_y_scale(cnt$n, axis = "x") +
    ggplot2::labs(title = "VEP consequence", x = "Variants", y = NULL) +
    theme_app(legend.position = "none")
  # Only override the default hue palette when an alternative is requested.
  if (!is.null(palette) && palette != "Default") {
    lv  <- levels(cnt$label)
    def <- stats::setNames(seq_along(lv), lv)
    p <- p + ggplot2::scale_fill_manual(values = apply_palette(def, palette))
  }
  p
}

# --- Top genes by sample count ----------------------------------------------
#' Horizontal bar chart of the genes carrying variants in the most samples.
#'   group_lookup: optional named character vector mapping family_id ->
#'                 diagnosis group ("MacTel", "HSAN1", "MacTel + HSAN1",
#'                 "Control"). When supplied AND two or more groups are present
#'                 in the data, bars are split into a stacked, diagnosis-coloured
#'                 chart; otherwise a single-colour bar is drawn (e.g. when only
#'                 MacTel patients are in view, there is nothing to distinguish).
plot_top_genes <- function(df, n_top = 25, group_lookup = NULL,
                           palette = "Default") {
  # Total distinct samples per gene: this picks and orders the top genes, and
  # (in the stacked case) is the total drawn at the end of each bar.
  totals <- df %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(n_samples = dplyr::n_distinct(family_id), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_samples)) %>%
    dplyr::slice_head(n = n_top)
  if (nrow(totals) == 0) return(NULL)
  gene_levels <- totals$SYMBOL[order(totals$n_samples)]   # ascending for the y axis

  # "Top N" only makes sense when the cap actually truncates the gene list; when
  # fewer genes than the cap are present, every gene is shown, so drop the count.
  gene_title <- if (nrow(totals) < n_top) "Genes by samples" else
    sprintf("Top %d genes by samples", n_top)

  # Which diagnosis groups are actually represented among these variants?
  groups_present <- character(0)
  if (!is.null(group_lookup)) {
    g <- unname(group_lookup[as.character(df$family_id)])
    groups_present <- sort(unique(g[!is.na(g)]))
  }

  # Single group (or no lookup): plain single-colour bar with count labels.
  if (length(groups_present) < 2) {
    d <- totals %>%
      dplyr::mutate(SYMBOL = factor(SYMBOL, levels = gene_levels))
    return(
      ggplot2::ggplot(d, ggplot2::aes(n_samples, SYMBOL)) +
        ggplot2::geom_col(fill = "#4C72B0", colour = "white", width = 0.7) +
        ggplot2::geom_text(ggplot2::aes(label = n_samples), hjust = -0.2,
                           size = 3, fontface = "bold") +
        ggplot2::scale_x_continuous(
          expand = ggplot2::expansion(mult = c(0, 0.12))) +
        ggplot2::labs(title = gene_title, x = "Samples", y = NULL) +
        theme_app(legend.position = "none",
                  axis.text.y = ggplot2::element_text(face = "italic"))
    )
  }

  # Two or more groups: stacked bar coloured by diagnosis. Count distinct
  # samples per (gene, group); each sample maps to exactly one group, so the
  # segments sum to the per-gene total drawn at the bar end.
  dd <- df %>%
    dplyr::filter(SYMBOL %in% totals$SYMBOL) %>%
    dplyr::mutate(diag_group = unname(group_lookup[as.character(family_id)])) %>%
    dplyr::filter(!is.na(diag_group)) %>%
    dplyr::distinct(SYMBOL, family_id, diag_group) %>%
    dplyr::count(SYMBOL, diag_group, name = "n_samples") %>%
    dplyr::mutate(
      SYMBOL     = factor(SYMBOL, levels = gene_levels),
      diag_group = factor(diag_group,
                          levels = intersect(names(COL_DIAG), groups_present)))
  totlab <- totals %>%
    dplyr::mutate(SYMBOL = factor(SYMBOL, levels = gene_levels))

  ggplot2::ggplot(dd, ggplot2::aes(n_samples, SYMBOL, fill = diag_group)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::geom_text(data = totlab, inherit.aes = FALSE,
                       ggplot2::aes(n_samples, SYMBOL, label = n_samples),
                       hjust = -0.2, size = 3, fontface = "bold") +
    ggplot2::scale_fill_manual(values = apply_palette(COL_DIAG, palette),
                               drop = TRUE, name = "Diagnosis") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(title = gene_title, x = "Samples", y = NULL) +
    theme_app(
      axis.text.y            = ggplot2::element_text(face = "italic"),
      legend.position        = "inside",
      legend.position.inside = c(0.98, 0.02),   # bottom-right, always clear here
      legend.justification   = c(1, 0),
      legend.background      = ggplot2::element_rect(
        fill = scales::alpha("white", 0.6), colour = NA))
}

# --- CADD vs REVEL scatter (interactive via plotly) --------------------------
# Clean, self-explanatory legend labels for the ClinVar colour channel.
CLNSIG_SCATTER_LABELS <- c(
  "Pathogenic"                   = "ClinVar P",
  "Pathogenic/Likely_pathogenic" = "ClinVar P/LP",
  "Likely_pathogenic"            = "ClinVar LP",
  "Conflicting_classifications"  = "ClinVar Conflicting",
  "Uncertain_significance"       = "ClinVar VUS",
  "Benign/Likely_benign"         = "ClinVar B/LB",
  "Not in ClinVar"               = "Not in ClinVar")

#' Colour encodes ClinVar classification (graded pathogenicity, blue when the
#' variant is not in ClinVar). Every point is a missense variant, so VEP impact
#' carries no information here and is not encoded. `key` carries the variant
#' identity so a plotly click can open its landing page.
plot_score_scatter <- function(df, threshold = 20) {
  d <- df %>%
    dplyr::filter(!is.na(CADD), !is.na(REVEL)) %>%
    dplyr::mutate(
      key = paste(CHROM, POS, REF, ALT),
      tooltip = sprintf("%s\n%s %s\nCADD %.1f | REVEL %.2f\n%s",
                        SYMBOL, HGVSc, ifelse(is.na(HGVSp_short), "", HGVSp_short),
                        CADD, REVEL, CLNSIG_clean))
  if (nrow(d) < 3) return(NULL)

  # Recode the colour variable to its display label rather than relabelling via
  # the scale's `labels=` argument: ggplotly ignores a manual scale's labels and
  # would fall back to the raw ClinVar terms in the interactive legend.
  raw_lvl  <- names(CLNSIG_SCATTER_LABELS)
  lbl      <- CLNSIG_SCATTER_LABELS[as.character(d$CLNSIG_clean)]
  lbl[is.na(lbl)] <- as.character(d$CLNSIG_clean)[is.na(lbl)]
  d$clin_label <- factor(lbl, levels = unique(CLNSIG_SCATTER_LABELS))

  # Graded ClinVar palette, but blue for unannotated variants (per request)
  # rather than the neutral grey used elsewhere. Key it by the display labels.
  clin_cols <- COL_CLNSIG
  clin_cols["Not in ClinVar"] <- "#1565C0"
  clin_cols <- stats::setNames(clin_cols[raw_lvl], CLNSIG_SCATTER_LABELS[raw_lvl])

  ggplot2::ggplot(d, ggplot2::aes(CADD, REVEL, colour = clin_label,
                                  text = tooltip, key = key)) +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_point(alpha = 0.8, size = 2.6) +
    ggplot2::scale_colour_manual(values = clin_cols,
                                 drop = TRUE, name = "ClinVar") +
    ggplot2::labs(title = "CADD vs REVEL (missense in silico)",
                  x = "CADD", y = "REVEL") +
    theme_app()
}

# --- Protein lollipop --------------------------------------------------------
#' Parse a 1-based amino-acid position out of an HGVSp string (e.g.
#' "p.Arg123Cys" -> 123). Returns NA when no residue number is present.
aa_position <- function(hgvsp) {
  suppressWarnings(as.integer(stringr::str_extract(hgvsp, "\\d+")))
}

#' Protein lollipop for one gene.
#'   gene_df  : all variant rows for the gene (needs HGVSp_short, CADD,
#'              CLNSIG_clean, CHROM, POS, REF, ALT, family_id)
#'   dom_df   : Pfam rows for the gene (SYMBOL, Protein_Length, Pfam, Domain,
#'              Start, End) from load_protein_domains(); may be empty/NULL
#'   gene     : gene symbol (for the title)
#'   sel_key  : "CHROM POS REF ALT" of the clicked variant to highlight (or NULL)
#'   label_all: if TRUE, label every lollipop point with its HGVSp (ggrepel).
#'              Used for the static gene report; the interactive view leaves it
#'              FALSE so the plot stays uncluttered.
#' Lollipop height = CADD, colour = ClinVar class, size = #samples carrying it.
#' Pfam domains are drawn as boxes on the protein backbone beneath the stems.
#'   italic_gene: if TRUE, italicise the gene symbol in the plot title via a
#'              plotmath expression. Left FALSE for the interactive view because
#'              ggplotly cannot convert plotmath titles (that view italicises the
#'              gene with an HTML tag after conversion instead).
plot_variant_lollipop <- function(gene_df, dom_df, gene, sel_key = NULL,
                                  label_all = FALSE, italic_gene = FALSE,
                                  threshold = 20, novel_keys = character(0),
                                  mark_novel = FALSE, affected_ids = NULL) {
  v <- gene_df %>%
    dplyr::mutate(
      aa  = aa_position(HGVSp_short),
      key = paste(CHROM, POS, REF, ALT)
    ) %>%
    dplyr::filter(!is.na(aa), !is.na(CADD))
  if (nrow(v) == 0) return(NULL)

  # Dot size counts only "affected" carriers (MacTel or HSAN1); controls still
  # appear in the total and on click but do not inflate the dot. When
  # affected_ids is not supplied, every carrier counts (backwards compatible).
  aff <- if (is.null(affected_ids)) unique(as.character(v$family_id))
         else as.character(affected_ids)

  # one lollipop per distinct variant
  vv <- v %>%
    dplyr::group_by(key, aa, CADD, CLNSIG_clean, HGVSp_short) %>%
    dplyr::summarise(
      n_carriers = dplyr::n_distinct(family_id),
      n_affected = dplyr::n_distinct(family_id[as.character(family_id) %in% aff]),
      .groups = "drop") %>%
    dplyr::mutate(tooltip = sprintf(
      "%s\nposition %d\nCADD %.1f\nClinVar: %s\ncarriers: %d (%d MacTel/HSAN1)",
      ifelse(is.na(HGVSp_short), "(no HGVSp)", HGVSp_short),
      aa, CADD, as.character(CLNSIG_clean), n_carriers, n_affected))

  # Variants flagged "novel for MacTel" are drawn as triangles when the toggle
  # is on. novel_keys use the CHROM||POS||REF||ALT form; vv$key is the
  # space-separated plotly click key, so match on the space form.
  vv$is_novel <- vv$key %in% gsub("||", " ", novel_keys, fixed = TRUE)
  show_novel   <- isTRUE(mark_novel) && any(vv$is_novel)

  # When no variant has more than one affected (MacTel/HSAN1) carrier the size
  # channel carries no information, so drop it (fixed dot size, no size legend).
  single_sample <- all(vv$n_affected <= 1)

  prot_len <- if (!is.null(dom_df) && nrow(dom_df) > 0)
    suppressWarnings(max(dom_df$Protein_Length, na.rm = TRUE)) else NA_real_
  if (!is.finite(prot_len)) prot_len <- max(vv$aa, na.rm = TRUE)

  ymax <- max(vv$CADD, na.rm = TRUE, 1)
  band <- ymax * 0.10                     # height of the backbone/domain band

  doms <- if (!is.null(dom_df)) dplyr::filter(dom_df, !is.na(Start), !is.na(End)) else dom_df[0, ]

  p <- ggplot2::ggplot() +
    # protein backbone
    ggplot2::annotate("segment", x = 1, xend = prot_len,
                      y = -band / 2, yend = -band / 2,
                      colour = "grey55", linewidth = 1.1)

  if (!is.null(doms) && nrow(doms) > 0) {
    p <- p +
      ggplot2::geom_rect(data = doms,
                         ggplot2::aes(xmin = Start, xmax = End,
                                      ymin = -band, ymax = 0, fill = Domain),
                         colour = "grey30", alpha = 0.9) +
      ggplot2::scale_fill_brewer(palette = "Set2", name = "Pfam domain")
  }

  # Point aesthetic: colour is always ClinVar; size is affected carriers unless
  # uninformative; shape splits known vs novel only when marking novel.
  point_aes <- if (single_sample) {
    if (show_novel)
      ggplot2::aes(x = aa, y = CADD, colour = CLNSIG_clean, shape = is_novel,
                   text = tooltip, key = key)
    else
      ggplot2::aes(x = aa, y = CADD, colour = CLNSIG_clean,
                   text = tooltip, key = key)
  } else {
    if (show_novel)
      ggplot2::aes(x = aa, y = CADD, colour = CLNSIG_clean, size = n_affected,
                   shape = is_novel, text = tooltip, key = key)
    else
      ggplot2::aes(x = aa, y = CADD, colour = CLNSIG_clean, size = n_affected,
                   text = tooltip, key = key)
  }

  p <- p +
    # stems + heads
    ggplot2::geom_segment(data = vv,
                          ggplot2::aes(x = aa, xend = aa, y = 0, yend = CADD),
                          colour = "grey70", linewidth = 0.5) +
    (if (single_sample)
       ggplot2::geom_point(data = vv, point_aes, size = 4)
     else
       ggplot2::geom_point(data = vv, point_aes)) +
    ggplot2::scale_colour_manual(values = COL_CLNSIG, drop = TRUE,
                                 name = "ClinVar") +
    (if (!single_sample)
       ggplot2::scale_size_continuous(
         range = c(2.5, 7), name = "Carriers",
         # integer breaks only — fractional carrier counts make no sense
         breaks = function(x) { b <- unique(round(scales::breaks_pretty(4)(x)))
                                b[b >= 1] })) +
    (if (show_novel)
       ggplot2::scale_shape_manual(
         values = c("FALSE" = 16, "TRUE" = 17),
         labels = c("FALSE" = "Known", "TRUE" = "Novel for MacTel"),
         name = "MacTel")) +
    ggplot2::geom_hline(yintercept = threshold, linetype = "dashed",
                        colour = "red", linewidth = 0.6)

  # highlight the clicked variant. If it has no amino-acid position it is not
  # protein-coding (e.g. intronic / splice / UTR) and cannot be drawn here, so
  # show a clear disclaimer instead — the variant detail above still applies.
  # The flag is returned as an attribute so the interactive (plotly) view can
  # add the same disclaimer (ggplot annotations are dropped by ggplotly).
  sel_not_coding <- FALSE
  if (!is.null(sel_key)) {
    sel <- dplyr::filter(vv, key == sel_key)
    if (nrow(sel) > 0) {
      p <- p +
        ggplot2::geom_point(data = sel,
                            ggplot2::aes(x = aa, y = CADD),
                            shape = 21, size = 8, stroke = 1.5,
                            colour = "black", fill = NA) +
        ggplot2::geom_text(data = sel,
                           ggplot2::aes(x = aa, y = CADD, label = HGVSp_short),
                           nudge_y = ymax * 0.09, vjust = 0,
                           fontface = "bold", size = 3.3)
    } else {
      sel_not_coding <- TRUE
      p <- p +
        ggplot2::annotate("label",
                          x = (1 + prot_len) / 2, y = ymax * 1.12,
                          label = "Selected variant is not protein coding",
                          fill = "#fff3cd", colour = "#664d03",
                          fontface = "bold", size = 3.5)
    }
  }

  # Label every point (gene report): HGVSp with leader lines, repelled so the
  # labels don't overlap even on gene-dense proteins.
  if (label_all && nrow(vv) > 0) {
    lab <- vv %>% dplyr::filter(!is.na(HGVSp_short), HGVSp_short != "")
    if (nrow(lab) > 0) {
      p <- p +
        ggrepel::geom_text_repel(
          data = lab,
          ggplot2::aes(x = aa, y = CADD, label = HGVSp_short),
          size = 2.6, fontface = "plain",
          min.segment.length = 0, max.overlaps = Inf,
          segment.size = 0.2, segment.colour = "grey60",
          box.padding = 0.3, point.padding = 0.1, seed = 1)
    }
  }

  p <- p +
    ggplot2::scale_x_continuous(limits = c(1, prot_len),
                                expand = ggplot2::expansion(mult = c(0.01, 0.03))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.24))) +
    ggplot2::labs(
      title    = if (isTRUE(italic_gene))
                   bquote(italic(.(gene)) * " protein lollipop")
                 else sprintf("%s protein lollipop", gene),
      subtitle = sprintf("%g aa | height = CADD (dashed = %g) | colour = ClinVar%s%s",
                         prot_len, threshold,
                         if (single_sample) "" else " | size = #MacTel/HSAN1 carriers",
                         if (show_novel) " | triangle = novel for MacTel" else ""),
      x = "Amino-acid position", y = "CADD") +
    theme_app()

  attr(p, "sel_not_coding") <- sel_not_coding
  p
}

# --- Pathway priority-variant summary ----------------------------------------
# A schematic figure: each metabolic pathway is a stacked lane, laid out
# left-to-right as a chain of gene nodes joined by reaction arrows. Every
# priority variant (>= 1 of: ClinVar P/LP, VEP HIGH, CADD >= threshold) is drawn
# as one marker above its gene. The three flags map to three independent visual
# channels so they can be read simultaneously:
#   * ClinVar P/LP  -> fill colour (crimson vs grey)
#   * VEP HIGH       -> shape       (triangle vs circle)
#   * CADD >= thr    -> border ring (black vs faint grey)
# Genes with no priority variant render as faded nodes, so absence is visible.
# The pathway/gene layout is data-driven from an editable spec (data/pathways.tsv)
# and the variant layer is recomputed from the data on every render.

plot_pathway_summary <- function(df, spec, threshold = 30) {
  if (is.null(spec) || nrow(spec) == 0) return(NULL)

  spec <- spec[!is.na(spec$symbol) & nzchar(spec$symbol), , drop = FALSE]
  spec$rank <- suppressWarnings(as.numeric(spec$rank))
  path_levels <- unique(spec$pathway)               # lane order = spec order
  n_lanes <- length(path_levels)
  # Top pathway gets the highest y so it sits at the top of the canvas.
  lane_y <- stats::setNames(rev(seq_len(n_lanes)) * 2, path_levels)
  spec <- spec[order(match(spec$pathway, path_levels), spec$rank), ]
  x_gap  <- 2.0
  spec$x <- spec$rank * x_gap
  spec$y <- lane_y[spec$pathway]

  # Legend level labels for the three independent flag channels.
  lv_clin <- c("ClinVar P/LP", "not P/LP")
  lv_high <- c("VEP HIGH", "not HIGH")
  lv_cadd <- c(sprintf("CADD >= %g", threshold), sprintf("CADD < %g", threshold))

  # --- priority-variant layer (one row per distinct variant) -----------------
  pv <- df
  if (!is.null(pv) && nrow(pv) > 0 && all(c("n_flags", "SYMBOL") %in% names(pv))) {
    pv <- pv[!is.na(pv$n_flags) & pv$n_flags >= 1, , drop = FALSE]
    pv <- dplyr::distinct(pv, SYMBOL, CHROM, POS, REF, ALT,
                          flag_clinvar, flag_high, flag_cadd)
    pv <- pv[pv$SYMBOL %in% spec$symbol, , drop = FALSE]
  } else {
    pv <- data.frame(SYMBOL = character(0), flag_clinvar = logical(0),
                     flag_high = logical(0), flag_cadd = logical(0))
  }
  pv$clin <- factor(ifelse(pv$flag_clinvar, lv_clin[1], lv_clin[2]), levels = lv_clin)
  pv$high <- factor(ifelse(pv$flag_high,    lv_high[1], lv_high[2]), levels = lv_high)
  pv$cadd <- factor(ifelse(pv$flag_cadd,    lv_cadd[1], lv_cadd[2]), levels = lv_cadd)

  # Stack each gene's variants in a compact centred grid just above its node.
  per_row <- 6; dx <- 0.19; dy <- 0.19; base_dy <- 0.44
  dots <- do.call(rbind, lapply(spec$symbol, function(sym) {
    v <- pv[pv$SYMBOL == sym, , drop = FALSE]
    if (nrow(v) == 0) return(NULL)
    node <- spec[spec$symbol == sym, ][1, ]
    m   <- nrow(v)
    idx <- seq_len(m) - 1L
    row <- idx %/% per_row
    xoff <- vapply(seq_len(m), function(i) {
      r <- row[i]; in_row <- which(row == r); pos <- match(i, in_row)
      (pos - 1 - (length(in_row) - 1) / 2) * dx
    }, numeric(1))
    data.frame(x = node$x + xoff, y = node$y + base_dy + row * dy,
               clin = v$clin, high = v$high, cadd = v$cadd,
               stringsAsFactors = FALSE)
  }))

  spec$n_pv <- vapply(spec$symbol, function(s) sum(pv$SYMBOL == s), integer(1))
  spec$has  <- spec$n_pv > 0

  # Intra-lane reaction arrows (consecutive gene nodes).
  seg <- do.call(rbind, lapply(path_levels, function(p) {
    s <- spec[spec$pathway == p, ]; s <- s[order(s$rank), ]
    if (nrow(s) < 2) return(NULL)
    data.frame(x = s$x[-nrow(s)], xend = s$x[-1], y = s$y[1], yend = s$y[1])
  }))

  x_max <- max(spec$x) + x_gap
  # Lane header labels sit at the far left, above each lane.
  lane_df <- data.frame(pathway = path_levels, y = lane_y[path_levels])

  p <- ggplot2::ggplot()
  # Faint lane bands.
  p <- p + ggplot2::geom_rect(
    data = lane_df,
    ggplot2::aes(xmin = 0.2, xmax = x_max, ymin = y - 0.8, ymax = y + 1.1),
    fill = "grey96", colour = NA)
  # Intra-lane reaction arrows.
  if (!is.null(seg)) p <- p + ggplot2::geom_segment(
    data = seg, ggplot2::aes(x = x + 0.42, xend = xend - 0.42, y = y, yend = yend),
    colour = "grey45", linewidth = 0.6,
    arrow = grid::arrow(length = grid::unit(6, "pt"), type = "closed"))
  # Gene nodes: filled/bold when they carry priority variants, faded when empty.
  # Drawn as fixed-colour layers (no colour aesthetic) so the colour scale is
  # free for the variant-dot border encoding below.
  spec_empty <- spec[!spec$has, , drop = FALSE]
  spec_has   <- spec[spec$has,  , drop = FALSE]
  if (nrow(spec_empty)) p <- p + ggplot2::geom_label(
    data = spec_empty, ggplot2::aes(x = x, y = y, label = label),
    fill = "white", colour = "grey65", fontface = "plain",
    linewidth = 0.3, label.r = grid::unit(4, "pt"), size = 3.3)
  if (nrow(spec_has)) p <- p + ggplot2::geom_label(
    data = spec_has, ggplot2::aes(x = x, y = y, label = label),
    fill = "white", colour = "grey10", fontface = "bold",
    linewidth = 0.5, label.r = grid::unit(4, "pt"), size = 3.5)
  # Priority-variant dots: three independent visual channels so each flag is
  # readable on its own — fill colour = ClinVar P/LP, shape = VEP HIGH impact,
  # border ring = CADD threshold.
  if (!is.null(dots) && nrow(dots) > 0) p <- p +
    ggplot2::geom_point(
      data = dots,
      ggplot2::aes(x = x, y = y, fill = clin, shape = high, colour = cadd),
      size = 3.4, stroke = 1.6) +
    # One legend block: show only the "included if" state of each channel; the
    # opposite (grey fill / circle / faint ring) is implicit. Two of the three
    # sub-guides carry no title so they read as a single labelled legend.
    ggplot2::scale_fill_manual(
      values = stats::setNames(c("#C62828", "#CFD8DC"), lv_clin),
      limits = lv_clin, breaks = lv_clin[1], drop = FALSE, name = "Included if") +
    ggplot2::scale_shape_manual(
      values = stats::setNames(c(24, 21), lv_high),
      limits = lv_high, breaks = lv_high[1], drop = FALSE, name = NULL) +
    ggplot2::scale_colour_manual(
      values = stats::setNames(c("#00BFA5", "grey75"), lv_cadd),
      limits = lv_cadd, breaks = lv_cadd[1], drop = FALSE, name = NULL) +
    ggplot2::guides(
      fill   = ggplot2::guide_legend(
        order = 1, override.aes = list(shape = 21, colour = "grey60")),
      shape  = ggplot2::guide_legend(
        order = 2, override.aes = list(fill = "grey75", colour = "grey60")),
      colour = ggplot2::guide_legend(
        order = 3, override.aes = list(shape = 21, fill = "grey85")))
  # Lane titles at the far left.
  p <- p + ggplot2::geom_text(
    data = lane_df,
    ggplot2::aes(x = 0.35, y = y + 1.35, label = pathway),
    hjust = 0, vjust = 1, fontface = "bold", size = 3.6, colour = "grey25")

  n_pv_total <- if (!is.null(dots)) nrow(dots) else 0
  p +
    ggplot2::coord_cartesian(clip = "off") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0.02, 0.05))) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0.08, 0.14))) +
    ggplot2::labs(
      title    = "Priority variants across serine / glycine / sphingolipid metabolism",
      subtitle = sprintf(
        "%d priority variant%s | grey / circle = criterion not met | faded gene = no priority variant",
        n_pv_total, if (n_pv_total == 1) "" else "s")) +
    theme_app(
      axis.title = ggplot2::element_blank(),
      axis.text  = ggplot2::element_blank(),
      axis.ticks = ggplot2::element_blank(),
      panel.grid = ggplot2::element_blank(),
      panel.border = ggplot2::element_blank(),
      legend.position = "right",
      legend.spacing.y = grid::unit(2, "pt"))
}

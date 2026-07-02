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

# --- IMPACT distribution -----------------------------------------------------
plot_impact <- function(df, palette = "Default") {
  df %>%
    dplyr::count(IMPACT, .drop = FALSE) %>%
    ggplot2::ggplot(ggplot2::aes(IMPACT, n, fill = IMPACT)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.4,
                       fontface = "bold", size = 3.5) +
    ggplot2::scale_fill_manual(values = apply_palette(COL_IMPACT, palette),
                               drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(title = "VEP impact", x = NULL, y = "Variants") +
    theme_app(legend.position = "none", axis.text.x = angle_x())
}

# --- Variant TYPE ------------------------------------------------------------
plot_type <- function(df, palette = "Default") {
  df %>%
    dplyr::count(TYPE, .drop = FALSE) %>%
    ggplot2::ggplot(ggplot2::aes(TYPE, n, fill = TYPE)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.4,
                       fontface = "bold", size = 3.5) +
    ggplot2::scale_fill_manual(values = apply_palette(COL_TYPE, palette),
                               drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.18))) +
    ggplot2::labs(title = "Variant type", x = NULL, y = "Variants") +
    theme_app(legend.position = "none", axis.text.x = angle_x())
}

# --- ClinVar -----------------------------------------------------------------
plot_clnsig <- function(df, palette = "Default") {
  cnt <- dplyr::count(df, CLNSIG_clean, .drop = FALSE)

  # Dynamic y-axis: when the counts span a wide range (e.g. Tier 2 genes where
  # almost everything is "Not in ClinVar" while the pathogenic categories have a
  # handful each), the small but informative bars vanish on a linear scale. If
  # the largest category is >= 20x the smallest non-zero one, switch to a
  # base-10 log scale so every category stays readable. A pseudo-log transform
  # is used so counts of 0 and 1 still render sensibly (plain log10 is undefined
  # at 0 and flattens 1 to the baseline). Otherwise use a linear axis.
  nz      <- cnt$n[cnt$n > 0]
  use_log <- length(nz) >= 2 && (max(nz) / min(nz)) >= 20
  expand  <- ggplot2::expansion(mult = c(0, 0.18))
  y_scale <- if (use_log) {
    ggplot2::scale_y_continuous(
      transform = scales::pseudo_log_trans(base = 10),
      breaks    = c(0, 1, 3, 10, 30, 100, 300, 1000),
      expand    = expand)
  } else {
    ggplot2::scale_y_continuous(expand = expand)
  }

  ggplot2::ggplot(cnt, ggplot2::aes(CLNSIG_clean, n, fill = CLNSIG_clean)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.4,
                       fontface = "bold", size = 3.2) +
    ggplot2::scale_fill_manual(values = apply_palette(COL_CLNSIG, palette),
                               drop = FALSE) +
    y_scale +
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
  ggplot2::ggplot(d, ggplot2::aes(CADD, fill = IMPACT)) +
    ggplot2::geom_histogram(binwidth = 2, colour = "white", alpha = 0.9) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed",
                        colour = "red", linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = apply_palette(COL_IMPACT, palette),
                               drop = TRUE, name = "Impact") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(title = "CADD distribution", x = "CADD", y = "Count") +
    theme_app()
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
  cnt$label <- factor(cnt$label, levels = cnt$label[order(cnt$n)])  # asc y axis

  p <- ggplot2::ggplot(cnt, ggplot2::aes(n, label, fill = label)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n), hjust = -0.2,
                       fontface = "bold", size = 3.5) +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
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
        ggplot2::labs(title = sprintf("Top %d genes by samples", n_top),
                      x = "Samples", y = NULL) +
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
    ggplot2::labs(title = sprintf("Top %d genes by samples", n_top),
                  x = "Samples", y = NULL) +
    theme_app(axis.text.y = ggplot2::element_text(face = "italic"))
}

# --- CADD vs REVEL scatter (interactive via plotly) --------------------------
plot_score_scatter <- function(df) {
  d <- df %>%
    dplyr::filter(!is.na(CADD), !is.na(REVEL)) %>%
    dplyr::mutate(
      tooltip = sprintf("%s\n%s %s\nCADD %.1f | REVEL %.2f\n%s",
                        SYMBOL, HGVSc, ifelse(is.na(HGVSp_short), "", HGVSp_short),
                        CADD, REVEL, CLNSIG_clean))
  if (nrow(d) < 3) return(NULL)
  ggplot2::ggplot(d, ggplot2::aes(CADD, REVEL, colour = IMPACT,
                                  shape = is_pathLP, text = tooltip)) +
    ggplot2::geom_point(alpha = 0.75, size = 2.4) +
    ggplot2::scale_colour_manual(values = COL_IMPACT, drop = TRUE,
                                 name = "Impact") +
    ggplot2::scale_shape_manual(values = c(`FALSE` = 16, `TRUE` = 17),
                                labels = c("Other", "ClinVar P/LP"),
                                name = "") +
    ggplot2::geom_hline(yintercept = 0.5, linetype = "dashed", colour = "grey60") +
    ggplot2::geom_vline(xintercept = 20, linetype = "dashed", colour = "grey60") +
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
                                  label_all = FALSE, italic_gene = FALSE) {
  v <- gene_df %>%
    dplyr::mutate(
      aa  = aa_position(HGVSp_short),
      key = paste(CHROM, POS, REF, ALT)
    ) %>%
    dplyr::filter(!is.na(aa), !is.na(CADD))
  if (nrow(v) == 0) return(NULL)

  # one lollipop per distinct variant; size by number of carriers
  vv <- v %>%
    dplyr::group_by(key, aa, CADD, CLNSIG_clean, HGVSp_short) %>%
    dplyr::summarise(n_carriers = dplyr::n_distinct(family_id),
                     .groups = "drop") %>%
    dplyr::mutate(tooltip = sprintf(
      "%s\nposition %d\nCADD %.1f\nClinVar: %s\nsamples: %d",
      ifelse(is.na(HGVSp_short), "(no HGVSp)", HGVSp_short),
      aa, CADD, as.character(CLNSIG_clean), n_carriers))

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

  p <- p +
    # stems + heads
    ggplot2::geom_segment(data = vv,
                          ggplot2::aes(x = aa, xend = aa, y = 0, yend = CADD),
                          colour = "grey70", linewidth = 0.5) +
    ggplot2::geom_point(data = vv,
                        ggplot2::aes(x = aa, y = CADD,
                                     colour = CLNSIG_clean, size = n_carriers,
                                     text = tooltip, key = key)) +
    ggplot2::scale_colour_manual(values = COL_CLNSIG, drop = TRUE,
                                 name = "ClinVar") +
    ggplot2::scale_size_continuous(range = c(2.5, 7), name = "Samples",
                                   breaks = scales::breaks_pretty(4)) +
    ggplot2::geom_hline(yintercept = 20, linetype = "dashed",
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
      subtitle = sprintf("%g aa | height = CADD (dashed = 20) | colour = ClinVar | size = #samples",
                         prot_len),
      x = "Amino-acid position", y = "CADD") +
    theme_app()

  attr(p, "sel_not_coding") <- sel_not_coding
  p
}

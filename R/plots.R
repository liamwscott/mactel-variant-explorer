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

# --- IMPACT distribution -----------------------------------------------------
plot_impact <- function(df) {
  df %>%
    dplyr::count(IMPACT, .drop = FALSE) %>%
    ggplot2::ggplot(ggplot2::aes(IMPACT, n, fill = IMPACT)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.4,
                       fontface = "bold", size = 3.5) +
    ggplot2::scale_fill_manual(values = COL_IMPACT, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(title = "VEP impact", x = NULL, y = "Variants") +
    theme_app(legend.position = "none")
}

# --- Variant TYPE ------------------------------------------------------------
plot_type <- function(df) {
  df %>%
    dplyr::count(TYPE, .drop = FALSE) %>%
    ggplot2::ggplot(ggplot2::aes(TYPE, n, fill = TYPE)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.4,
                       fontface = "bold", size = 3.5) +
    ggplot2::scale_fill_manual(values = COL_TYPE, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(title = "Variant type", x = NULL, y = "Variants") +
    theme_app(legend.position = "none")
}

# --- ClinVar -----------------------------------------------------------------
plot_clnsig <- function(df) {
  df %>%
    dplyr::count(CLNSIG_clean, .drop = FALSE) %>%
    ggplot2::ggplot(ggplot2::aes(CLNSIG_clean, n, fill = CLNSIG_clean)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n), vjust = -0.4,
                       fontface = "bold", size = 3.2) +
    ggplot2::scale_fill_manual(values = COL_CLNSIG, drop = FALSE) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::scale_x_discrete(labels = function(x)
      stringr::str_wrap(gsub("_", " ", x), 12)) +
    ggplot2::labs(title = "ClinVar classification", x = NULL, y = "Variants") +
    theme_app(legend.position = "none",
              axis.text.x = ggplot2::element_text(size = 8))
}

# --- CADD histogram ----------------------------------------------------------
plot_cadd <- function(df, threshold = 20) {
  d <- dplyr::filter(df, !is.na(CADD))
  if (nrow(d) == 0) return(NULL)
  ggplot2::ggplot(d, ggplot2::aes(CADD, fill = IMPACT)) +
    ggplot2::geom_histogram(binwidth = 2, colour = "white", alpha = 0.9) +
    ggplot2::geom_vline(xintercept = threshold, linetype = "dashed",
                        colour = "red", linewidth = 0.8) +
    ggplot2::scale_fill_manual(values = COL_IMPACT, drop = FALSE,
                               name = "Impact") +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.08))) +
    ggplot2::labs(title = "CADD distribution",
                  subtitle = sprintf("dashed line = CADD %g", threshold),
                  x = "CADD", y = "Count") +
    theme_app()
}

# --- Inheritance -------------------------------------------------------------
plot_inheritance <- function(df) {
  df %>%
    dplyr::count(inheritance) %>%
    ggplot2::ggplot(ggplot2::aes(stats::reorder(inheritance, n), n,
                                 fill = inheritance)) +
    ggplot2::geom_col(colour = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n), hjust = -0.2,
                       fontface = "bold", size = 3.5) +
    ggplot2::scale_y_continuous(expand = ggplot2::expansion(mult = c(0, 0.15))) +
    ggplot2::coord_flip() +
    ggplot2::labs(title = "Inheritance mode", x = NULL, y = "Variants") +
    theme_app(legend.position = "none")
}

# --- Top genes by family count ----------------------------------------------
plot_top_genes <- function(df, n_top = 25) {
  d <- df %>%
    dplyr::group_by(SYMBOL) %>%
    dplyr::summarise(n_samples  = dplyr::n_distinct(family_id),
                     n_variants = dplyr::n(), .groups = "drop") %>%
    dplyr::arrange(dplyr::desc(n_samples)) %>%
    dplyr::slice_head(n = n_top)
  if (nrow(d) == 0) return(NULL)
  d %>%
    dplyr::mutate(SYMBOL = forcats::fct_reorder(SYMBOL, n_samples)) %>%
    ggplot2::ggplot(ggplot2::aes(n_samples, SYMBOL)) +
    ggplot2::geom_col(fill = "#4C72B0", colour = "white", width = 0.7) +
    ggplot2::geom_text(ggplot2::aes(label = n_samples), hjust = -0.2,
                       size = 3, fontface = "bold") +
    ggplot2::scale_x_continuous(expand = ggplot2::expansion(mult = c(0, 0.12))) +
    ggplot2::labs(title = sprintf("Top %d genes by samples", n_top),
                  x = "Samples", y = NULL) +
    theme_app()
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
    ggplot2::scale_colour_manual(values = COL_IMPACT, drop = FALSE,
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

# =============================================================================
# R/load_data.R — load & clean a Cavalier candidate-variant CSV
# Expected columns: the 71-column Cavalier output (family_id, CHROM, POS, REF,
# ALT, CADD, CLNSIG, IMPACT, TYPE, SYMBOL, HGVSc, HGVSp, am_class,
# am_pathogenicity, REVEL, SpliceAI_max, gnomad_AF, inheritance, variant_id, …)
# =============================================================================

# Canonical orderings / palettes reused across the app -----------------------
IMPACT_LEVELS <- c("HIGH", "MODERATE", "LOW", "MODIFIER")
TYPE_LEVELS   <- c("LOF", "SPLICING", "MISSENSE", "OTHER")

CLNSIG_LEVELS <- c("Pathogenic", "Likely_pathogenic", "Pathogenic/Likely_pathogenic",
                   "Conflicting_classifications", "Uncertain_significance",
                   "Benign/Likely_benign", "Not in ClinVar")

COL_IMPACT <- c(HIGH = "#D62728", MODERATE = "#FF7F0E",
                LOW = "#2CA02C", MODIFIER = "#AEC7E8")

COL_TYPE <- c(LOF = "#9467BD", SPLICING = "#E377C2",
              MISSENSE = "#FF7F0E", OTHER = "#BCBD22")

COL_CLNSIG <- c("Pathogenic"                   = "#D62728",
                "Likely_pathogenic"             = "#FF7F0E",
                "Pathogenic/Likely_pathogenic"  = "#9467BD",
                "Conflicting_classifications"   = "#BCBD22",
                "Uncertain_significance"        = "#AEC7E8",
                "Benign/Likely_benign"          = "#2CA02C",
                "Not in ClinVar"                = "#DDDDDD")

#' Clean a raw Cavalier dataframe into the structure the app expects.
clean_variants <- function(df) {
  stopifnot(is.data.frame(df))

  required <- c("family_id", "SYMBOL", "CHROM", "POS", "REF", "ALT",
                "CADD", "CLNSIG", "IMPACT", "TYPE")
  missing  <- setdiff(required, names(df))
  if (length(missing) > 0) {
    stop("Input is missing required columns: ", paste(missing, collapse = ", "))
  }

  num_cols <- intersect(c("CADD", "REVEL", "am_pathogenicity", "SpliceAI_max",
                          "gnomad_AF", "gnomad_AC", "phyloP100", "POS",
                          "AF", "AC", "AN", "SIFT_score", "PolyPhen_score"),
                        names(df))

  df <- df %>%
    dplyr::mutate(dplyr::across(dplyr::all_of(num_cols),
                                ~ suppressWarnings(as.numeric(.)))) %>%
    dplyr::mutate(
      CLNSIG_clean = dplyr::case_when(
        grepl("Pathogenic/Likely_pathogenic", CLNSIG) ~ "Pathogenic/Likely_pathogenic",
        grepl("^Pathogenic$", CLNSIG)                 ~ "Pathogenic",
        grepl("^Likely_pathogenic$", CLNSIG)          ~ "Likely_pathogenic",
        grepl("Uncertain", CLNSIG)                    ~ "Uncertain_significance",
        grepl("Conflicting", CLNSIG)                  ~ "Conflicting_classifications",
        grepl("[Bb]enign", CLNSIG)                    ~ "Benign/Likely_benign",
        is.na(CLNSIG) | CLNSIG %in% c("NA", "")       ~ "Not in ClinVar",
        TRUE                                          ~ CLNSIG
      ),
      IMPACT       = factor(IMPACT,       levels = IMPACT_LEVELS),
      TYPE         = factor(TYPE,         levels = TYPE_LEVELS),
      CLNSIG_clean = factor(CLNSIG_clean, levels = CLNSIG_LEVELS),
      inheritance  = ifelse(is.na(inheritance) | inheritance == "",
                            "unknown", inheritance),
      HGVSp_short  = stringr::str_extract(HGVSp, "p\\..*"),
      is_pathLP    = CLNSIG_clean %in% c("Pathogenic", "Likely_pathogenic",
                                         "Pathogenic/Likely_pathogenic")
    )

  df
}

#' Load and clean from a file path.
load_variants <- function(path) {
  df <- readr::read_csv(path, show_col_types = FALSE, guess_max = 5000)
  clean_variants(df)
}

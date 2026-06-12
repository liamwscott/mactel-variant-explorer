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
      # Strip the trailing "RLA" suffix from sample IDs for display
      family_id    = stringr::str_remove(as.character(family_id), "RLA$"),
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

# Known monogenic MacTel genes — fallback Tier 1 if no lookup file is present.
TIER1_FALLBACK <- c("PHGDH", "SPTLC1", "SPTLC2")

#' Load a gene -> tier lookup (TSV with columns Gene_Symbol, Tier).
#' Returns a tibble with SYMBOL + Tier ("Tier 1"/"Tier 2"), or NULL if absent.
load_gene_tiers <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  t <- readr::read_tsv(path, show_col_types = FALSE)
  names(t)[1:2] <- c("SYMBOL", "Tier")
  t %>%
    dplyr::mutate(
      Tier = paste0("Tier ", gsub("[^0-9]", "", as.character(Tier)))
    ) %>%
    dplyr::select(SYMBOL, Tier) %>%
    dplyr::distinct()
}

#' Load the gene-information table (TSV with columns Gene_Symbol, Tier,
#' Ensembl_ID, Chromosome, Evidence_Category, Evidence_Detail, Gene_Description).
#' Returns a tibble keyed by SYMBOL, or NULL if the file is absent.
load_gene_info <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  g <- readr::read_tsv(path, show_col_types = FALSE)
  names(g)[names(g) == "Gene_Symbol"] <- "SYMBOL"
  g %>%
    dplyr::mutate(
      Tier = paste0("Tier ", gsub("[^0-9]", "", as.character(Tier)))
    ) %>%
    dplyr::distinct(SYMBOL, .keep_all = TRUE)
}

#' Load the bundled Pfam protein-domain table (TSV produced offline by
#' scripts/fetch_protein_domains.py). Columns: Gene_Symbol, UniProt,
#' Protein_Length, Pfam, Domain, Start, End. Returns a tibble keyed by SYMBOL
#' with numeric Start/End/Protein_Length, or NULL if the file is absent.
load_protein_domains <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  d <- readr::read_tsv(path, show_col_types = FALSE)
  names(d)[names(d) == "Gene_Symbol"] <- "SYMBOL"
  d %>%
    dplyr::mutate(
      Protein_Length = suppressWarnings(as.numeric(Protein_Length)),
      Start          = suppressWarnings(as.numeric(Start)),
      End            = suppressWarnings(as.numeric(End))
    )
}

# Data-group flag columns (value 1 = member) -> display label for the
# sample-explorer tags. Mito_haplo is deliberately excluded.
SAMPLE_TAG_COLS <- c(
  WES            = "WES",
  Golden_cohort  = "Golden cohort",
  Clinical_trial = "Clinical trial",
  HSAN1_variant  = "HSAN1",
  Early_onset    = "Early onset",
  Low_PRS        = "Low PRS",
  Chr_5          = "Chr 5",
  Family_1       = "Family 1",
  Other          = "Other"
)

#' Load the per-sample information sheet (TSV).
#' Accepts either the real sheet (keyed by `Manifest_Sample_ID`, e.g. A0001RLA)
#' or a de-identified sheet (already keyed by `family_id`, e.g. FAMILY001).
#' Returns a tibble keyed by `family_id` (matching the variant data) with logical
#' case/control flags plus the raw group flags used for tags, or NULL if absent.
load_sample_info <- function(path) {
  if (is.null(path) || !file.exists(path)) return(NULL)
  s <- readr::read_tsv(path, show_col_types = FALSE)
  if (!("family_id" %in% names(s)) && "Manifest_Sample_ID" %in% names(s)) {
    s$family_id <- stringr::str_remove(as.character(s$Manifest_Sample_ID), "RLA$")
  }
  s %>%
    dplyr::mutate(
      family_id  = as.character(family_id),
      is_mactel  = !is.na(MacTel_Diagnosis) & MacTel_Diagnosis == "yes",
      is_hsan1   = !is.na(HSAN1_variant) & as.numeric(HSAN1_variant) == 1,
      is_control = !is_mactel & !is_hsan1
    )
}

#' Add a Tier column to a variant dataframe.
#' Uses tier_df if supplied; otherwise falls back to the hardcoded Tier 1 list.
#' Genes absent from the lookup become "Unassigned".
annotate_tier <- function(df, tier_df = NULL) {
  if (is.null(tier_df)) {
    df$Tier <- ifelse(df$SYMBOL %in% TIER1_FALLBACK, "Tier 1", "Tier 2")
    return(df)
  }
  df %>%
    dplyr::left_join(tier_df, by = "SYMBOL") %>%
    dplyr::mutate(Tier = dplyr::if_else(is.na(Tier), "Unassigned", Tier))
}

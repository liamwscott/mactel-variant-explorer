# =============================================================================
# R/alphafold.R — on-demand AlphaFold structure fetching for the 3D viewer
#
# Structures are pulled from the AlphaFold DB the first time a gene is viewed
# and cached on disk, so the app stays light (nothing bundled in the repo) and
# subsequent views are instant. The 3D viewer needs internet the first time a
# given structure is shown; every other feature works fully offline.
# =============================================================================

# Per-residue pLDDT confidence colours (AlphaFold convention): the model stores
# pLDDT in the B-factor field, so the viewer colours by atom.b.
#   >90 very high (dark blue) · 70-90 confident (cyan) ·
#   50-70 low (yellow) · <50 very low / likely disordered (orange)
AF_PLDDT_COLORFUNC <- paste0(
  "function(atom){var b=atom.b;",
  "return b>90?'#0053D6':b>70?'#65CBF3':b>50?'#FFDB13':'#FF7D45';}")

#' Writable on-disk cache directory for downloaded AlphaFold models.
#' Uses the per-user R cache dir (R >= 4.0); falls back to the session tempdir
#' when that is unavailable (e.g. read-only install).
af_cache_dir <- function() {
  d <- tryCatch(tools::R_user_dir("MacTelVariantExplorer", "cache"),
                error = function(e) file.path(tempdir(), "mactel_cache"))
  d <- file.path(d, "alphafold")
  if (!dir.exists(d)) dir.create(d, recursive = TRUE, showWarnings = FALSE)
  d
}

#' Resolve a gene symbol to a UniProt accession using the bundled protein-domain
#' table (which already carries one reviewed human accession per gene).
#' Returns a single accession string, or NA_character_ when none is mapped.
af_uniprot_for_gene <- function(symbol, domains) {
  if (is.null(domains) || is.null(symbol) || length(symbol) != 1 || is.na(symbol))
    return(NA_character_)
  if (!all(c("SYMBOL", "UniProt") %in% names(domains))) return(NA_character_)
  u <- domains$UniProt[domains$SYMBOL == symbol]
  u <- u[!is.na(u) & nzchar(u)]
  if (length(u)) u[1] else NA_character_
}

#' Ask the AlphaFold DB API for the canonical PDB URL of an accession. Going
#' through the API (rather than guessing a model version) keeps the app working
#' as the DB bumps versions over time. PDB is used rather than mmCIF because
#' 3Dmol reliably reads per-residue pLDDT from the PDB B-factor column, whereas
#' the mmCIF B_iso field is not surfaced as atom.b. Returns a URL or NULL.
af_pdb_url <- function(uniprot) {
  url <- sprintf("https://alphafold.ebi.ac.uk/api/prediction/%s", uniprot)
  tryCatch({
    resp <- httr::GET(url, httr::timeout(25))
    if (httr::status_code(resp) != 200) return(NULL)
    meta <- httr::content(resp, as = "parsed", type = "application/json")
    if (length(meta) == 0) return(NULL)
    u <- meta[[1]]$pdbUrl
    if (is.null(u) || !nzchar(u)) NULL else u
  }, error = function(e) NULL)
}

#' Fetch the AlphaFold model for a UniProt accession as PDB text, caching it on
#' disk (keyed by accession). The PDB B-factor column carries per-residue pLDDT.
#' Returns a single string of PDB content, or NULL when the accession is missing,
#' the protein is not modelled, or the download fails (offline). The cached copy
#' is reused on every subsequent view.
fetch_alphafold_pdb <- function(uniprot) {
  if (is.null(uniprot) || length(uniprot) != 1 ||
      is.na(uniprot) || !nzchar(uniprot)) return(NULL)

  cache <- file.path(af_cache_dir(), sprintf("AF-%s-F1.pdb", uniprot))
  if (file.exists(cache) && file.info(cache)$size > 0) {
    return(paste(readLines(cache, warn = FALSE), collapse = "\n"))
  }

  url <- af_pdb_url(uniprot)
  if (is.null(url)) return(NULL)
  tryCatch({
    resp <- httr::GET(url, httr::timeout(25))
    if (httr::status_code(resp) != 200) return(NULL)
    txt <- httr::content(resp, as = "text", encoding = "UTF-8")
    if (is.null(txt) || !nzchar(txt) ||
        !(startsWith(txt, "HEADER") || grepl("\nATOM ", txt, fixed = TRUE)))
      return(NULL)
    writeLines(txt, cache)        # cache for next time
    txt
  }, error = function(e) NULL)
}

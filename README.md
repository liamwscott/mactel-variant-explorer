# MacTel Variant Explorer

An interactive **R Shiny** app for exploring Cavalier candidate-variant output
(rare-variant prioritisation in the MacTel Tier 1/2 gene list).

Point it at a Cavalier CSV and interactively filter, plot, and export variants
by gene, VEP impact, variant type, ClinVar class, inheritance mode, CADD,
REVEL, gnomAD frequency, and combined "priority flags".

---

## Quick start

```r
# from R / RStudio
shiny::runApp("path/to/mactel_variant_explorer")
```

or from the shell:

```bash
Rscript -e 'shiny::runApp("path/to/mactel_variant_explorer", launch.browser = TRUE)'
```

The app opens in your browser. On startup it loads, in order of preference:

1. `data/candidate_variants.csv` — your real data (git-ignored, never committed)
2. `data/example_variants.csv` — de-identified example shipped with the repo
3. any CSV you upload via the sidebar **Data** panel

---

## Required packages

```r
install.packages(c(
  "shiny", "bslib", "bsicons", "DT", "plotly",
  "ggplot2", "dplyr", "tidyr", "readr", "stringr",
  "forcats", "scales"
))
```

---

## Input format

Expects the **71-column Cavalier output** structure. The app strictly requires:

`family_id, SYMBOL, CHROM, POS, REF, ALT, CADD, CLNSIG, IMPACT, TYPE`

and makes use of (when present):
`HGVSc, HGVSp, am_class, am_pathogenicity, REVEL, SpliceAI_max, gnomad_AF,
inheritance, variant_id`.

---

## Features

| Tab | What it does |
|-----|--------------|
| **Overview** | Live summary plots: VEP impact, variant type, inheritance, CADD distribution, ClinVar classification, top genes by family count |
| **Variant table** | Searchable/sortable `DT` table of all filtered variants, ClinVar & HIGH-impact rows highlighted, CSV export |
| **Score scatter** | Interactive `plotly` CADD vs REVEL plot; hover shows gene + HGVS + ClinVar |
| **Priority variants** | Variants meeting ≥ N of {ClinVar P/LP, HIGH impact, CADD ≥ threshold}, with a plain-English "why prioritised" column |
| **Gene summary** | Per-gene rollup: variants, families, P/LP count, HIGH count, max CADD/REVEL |

All sidebar filters apply globally across every tab. Value boxes at the top show
live variant / gene / family / ClinVar-P-LP counts for the current filter set.

---

## Privacy

⚠️ **`family_id` contains real family identifiers.**
`data/candidate_variants.csv` and any `*.real.csv` are **git-ignored** and must
never be committed. Only the de-identified `data/example_variants.csv`
(family IDs replaced with `FAMILY001…`) lives in the repo.

---

## Layout

```
mactel_variant_explorer/
├── app.R                       # UI + server
├── R/
│   ├── load_data.R             # CSV load + cleaning, palettes, level orders
│   └── plots.R                 # reusable ggplot builders
├── data/
│   ├── candidate_variants.csv  # REAL data (git-ignored)
│   └── example_variants.csv    # de-identified example (committed)
├── .gitignore
└── README.md
```

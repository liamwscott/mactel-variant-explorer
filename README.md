# MacTel Variant Explorer

An interactive **R Shiny** app for exploring Cavalier candidate-variant output
(rare-variant prioritisation in the MacTel Tier 1/2 gene list).

Point it at a Cavalier CSV and interactively filter, plot, and export variants
by gene, VEP impact, variant type, ClinVar class, inheritance mode, CADD,
REVEL, gnomAD frequency, and combined "priority flags".

---

## Quick start

1. **Install R** (once): https://cran.r-project.org
2. **Unzip** the `mactel_variant_explorer` folder.
3. **Double-click the launcher** for your system:
   - **macOS** → `MacTel Explorer (Mac).app`
   - **Windows** → `Run MacTel Explorer (Windows).bat`

The app opens in your browser. The first launch installs the required R packages
automatically; later launches start in seconds. To stop the app, close the
terminal window the launcher opened.

> **macOS first-run note:** macOS blocks downloaded launchers it can't verify
> (you'll see *"…command" Not Opened*). **Do not click "Move to Trash."** Instead:
>
> 1. Click **Done**.
> 2. Open **System Settings → Privacy & Security**, scroll to the **Security**
>    section, and click **Open Anyway** next to the blocked launcher.
> 3. Authenticate, then click **Open Anyway** again. It won't ask after this.
>
> Or, in Terminal, clear the quarantine flag on the unzipped folder once:
> `xattr -dr com.apple.quarantine "/path/to/mactel_variant_explorer"`

---

## Running in RStudio

Open the project folder in RStudio and run:

```r
shiny::runApp("path/to/mactel_variant_explorer/app")
```

On startup the app loads, in order of preference:

1. `app/data/candidate_variants.csv` — your real data (git-ignored, never committed)
2. any CSV you upload via the sidebar **Data** panel

No participant data ships with the repo; supply your own Cavalier CSV.

---

## Required packages

```r
install.packages(c(
  "shiny", "bslib", "bsicons", "DT", "plotly",
  "ggplot2", "dplyr", "tidyr", "readr", "stringr",
  "forcats", "scales", "shinyFiles", "jsonlite",
  "r3dmol", "httr", "patchwork"
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

## Layout

```
mactel_variant_explorer/
├── MacTel Explorer (Mac).app       # macOS double-click launcher (DNA app icon)
├── Run MacTel Explorer (Windows).bat   # Windows double-click launcher
├── README.md
├── .gitignore
└── app/                            # all the app's code + data lives here
    ├── Run MacTel Explorer (Mac).command  # script the .app runs (also works on its own)
    ├── launch.R                    # bootstrap: installs deps + starts the app
    ├── app.R                       # UI + server
    ├── R/
    │   ├── load_data.R             # CSV load + cleaning, palettes, level orders
    │   └── plots.R                 # reusable ggplot builders
    ├── scripts/
    │   └── fetch_protein_domains.py  # dev tool: regenerate protein_domains.tsv
    └── data/
        └── candidate_variants.csv  # REAL data (git-ignored, supplied separately)
```

#!/usr/bin/env bash
#
# build_igv_reports_folder.sh
# -----------------------------------------------------------------------------
# Build a slimmed, shippable IGV-reports folder for the MacTel Variant Explorer.
#
# A full per-sample output folder (one sub-folder per individual) also contains
# large .pptx / .pdf / .png / .csv.gz files the app never uses. This script
# copies ONLY the "<SAMPLE>.igv_report.html" files, preserving the per-sample
# sub-folder layout the app expects:
#
#     <dest>/<SAMPLE>/<SAMPLE>.igv_report.html
#
# Point the app's "IGV reports folder" picker at <dest>.
#
# Usage:
#   tools/build_igv_reports_folder.sh <source_folder> [dest_folder]
#
# Examples:
#   tools/build_igv_reports_folder.sh ~/Downloads/by_family
#   tools/build_igv_reports_folder.sh ~/Downloads/by_family ./igv_reports_slim
#
# Notes:
#   - Reports render via an online igv.js CDN, so collaborators need internet.
#   - Only samples that actually have a report are copied (others are skipped).
# -----------------------------------------------------------------------------
set -euo pipefail

SRC="${1:-}"
DEST="${2:-./igv_reports}"

if [[ -z "$SRC" || ! -d "$SRC" ]]; then
  echo "Usage: $0 <source_folder> [dest_folder]" >&2
  echo "  Copies only *.igv_report.html (with their per-sample sub-folders)" >&2
  echo "  into <dest_folder> (default: ./igv_reports) for sharing." >&2
  exit 1
fi

mkdir -p "$DEST"

count=0
while IFS= read -r -d '' f; do
  sample="$(basename "$(dirname "$f")")"
  mkdir -p "$DEST/$sample"
  cp "$f" "$DEST/$sample/"
  count=$((count + 1))
done < <(find "$SRC" -name "*.igv_report.html" -print0)

echo "Copied ${count} IGV report(s) into: ${DEST}"
if command -v du >/dev/null 2>&1; then
  echo "Total size: $(du -sh "$DEST" | cut -f1)"
fi

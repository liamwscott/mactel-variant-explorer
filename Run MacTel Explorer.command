#!/usr/bin/env bash
# =============================================================================
# Run MacTel Explorer.command — macOS double-click launcher.
#
# Double-click this file in Finder to start the MacTel Variant Explorer.
# (The first time, macOS may ask you to confirm running it — see SETUP below.)
# =============================================================================

# Move into the folder this script lives in, so R finds app.R and data/.
DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$DIR" || exit 1

# Locate Rscript: PATH first, then the usual install locations.
if command -v Rscript >/dev/null 2>&1; then
  RSCRIPT="$(command -v Rscript)"
elif [ -x /usr/local/bin/Rscript ]; then
  RSCRIPT=/usr/local/bin/Rscript
elif [ -x /opt/homebrew/bin/Rscript ]; then
  RSCRIPT=/opt/homebrew/bin/Rscript
elif [ -x /Library/Frameworks/R.framework/Resources/bin/Rscript ]; then
  RSCRIPT=/Library/Frameworks/R.framework/Resources/bin/Rscript
else
  osascript -e 'display alert "R is not installed" message "Please install R from https://cran.r-project.org first, then double-click this launcher again."'
  echo "R not found. Install it from https://cran.r-project.org and try again."
  read -r -p "Press Return to close..."
  exit 1
fi

echo "Using R at: $RSCRIPT"
"$RSCRIPT" launch.R

# Keep the window open if something went wrong so the user can read the error.
status=$?
if [ $status -ne 0 ]; then
  echo ""
  echo "The app exited with an error (code $status)."
  read -r -p "Press Return to close..."
fi

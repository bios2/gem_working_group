#!/usr/bin/env bash
# init.sh
# Run once after cloning this repository.
# Clones external dependencies that are not tracked in git.

set -euo pipefail

REPO_ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONTEXT_DIR="$REPO_ROOT/experiments/sunday_vegetation_experiments_design/context"
MADINGLEY_DIR="$CONTEXT_DIR/MadingleyR"
MADINGLEY_URL="https://github.com/MadingleyR/MadingleyR.git"

echo "GEM Working Group -- environment init"
echo "======================================"

# Clone MadingleyR if not already present
if [ -d "$MADINGLEY_DIR/.git" ]; then
    echo "MadingleyR already cloned at:"
    echo "  $MADINGLEY_DIR"
    echo "Pulling latest changes..."
    git -C "$MADINGLEY_DIR" pull --ff-only
else
    echo "Cloning MadingleyR..."
    mkdir -p "$CONTEXT_DIR"
    git clone "$MADINGLEY_URL" "$MADINGLEY_DIR"
    echo "Done."
fi

echo ""
echo "MadingleyR version: $(git -C "$MADINGLEY_DIR" describe --tags --always)"
echo ""
echo "Next steps:"
echo "  Install the R package (run in R):"
echo "    remotes::install_local('experiments/sunday_vegetation_experiments_design/context/MadingleyR/Package')"
echo ""
echo "Init complete."

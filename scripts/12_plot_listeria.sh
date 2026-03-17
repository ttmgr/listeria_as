#!/bin/bash
# -----------------------------------------------------------------------------
# Optional plot step: make quick Listeria read-count and base-count figures.
# Input: processing/listeria/listeria_summary.tsv
# Output: processing/listeria/plots/*.png and *.pdf
# Run: sbatch scripts/12_plot_listeria.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

INPUT="${WORK_DIR}/processing/listeria/listeria_summary.tsv"
OUTPUT="${WORK_DIR}/processing/listeria/plots"
echo "Plotting Listeria results..."
echo "Start time: $(date)"
python3 "${SCRIPT_DIR}/plot_listeria.py" "$INPUT" "$OUTPUT"
echo "Finished: $(date)"

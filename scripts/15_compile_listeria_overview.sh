#!/bin/bash
# -----------------------------------------------------------------------------
# Step 15: Build one integrated Listeria overview table and plots.
# Input: Listeria read summary, Listeria contig summary, and NanoStat outputs
# Output: processing/listeria/overview/listeria_overview.csv (+ figures)
# Run: sbatch --dependency=afterok:<LISTERIA_READS>:<LISTERIA_CONTIGS>:<NANOSTAT> scripts/15_compile_listeria_overview.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

echo "Compiling comprehensive Listeria overview..."
echo "Start time: $(date)"
python3 "${SCRIPT_DIR}/15_compile_listeria_overview.py" "$WORK_DIR"
echo "Finished: $(date)"

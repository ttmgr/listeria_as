#!/bin/bash
# -----------------------------------------------------------------------------
# Step 16: Merge AMRFinder outputs into summary tables for reporting.
# Input: processing/amrfinder/{reads,flye,mdbg,myloasm}/amrfinder_*.tsv
# Output: processing/amrfinder/overview/*.csv
# Run: sbatch --dependency=afterok:<AMRFINDER_JOB> scripts/16_compile_amr_overview.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

echo "Compiling AMRFinderPlus overview..."
echo "Start time: $(date)"
python3 "${SCRIPT_DIR}/16_compile_amr_overview.py" "$WORK_DIR"
echo "Finished: $(date)"

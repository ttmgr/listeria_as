#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# -----------------------------------------------------------------------------
# Step 16: Merge AMRFinder outputs into summary tables for reporting.
# Input: processing/amrfinder/{reads,flye,mdbg,myloasm}/amrfinder_*.tsv
# Output: processing/amrfinder/overview/*.csv
# Run: sbatch --dependency=afterok:<AMRFINDER_JOB> scripts/16_compile_amr_overview.sh
# -----------------------------------------------------------------------------
WORK_DIR="/path/to/project"
echo "Compiling AMRFinderPlus overview..."
echo "Start time: $(date)"
python3 "${WORK_DIR}/scripts/16_compile_amr_overview.py" "$WORK_DIR"
echo "Finished: $(date)"

#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# -----------------------------------------------------------------------------
# Step 19: Build a fast reads-focused report before assembly-heavy steps finish.
# Input: read metrics, Listeria read summary, and AMR read summary
# Output: processing/report/reads_report.html
# Run: sbatch scripts/19_reads_report.sh
# -----------------------------------------------------------------------------
WORK_DIR="/path/to/project"
echo "Generating reads-only report..."
echo "Start time: $(date)"
python3 "${WORK_DIR}/scripts/19_reads_report.py" "$WORK_DIR"
echo "Finished: $(date)"
echo "Report: ${WORK_DIR}/processing/report/reads_report.html"

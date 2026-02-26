#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# ============================================================
# Step 19: Reads-Only QC Report (Reads + Listeria + AMR)
# Submit as: sbatch 19_reads_report.sh
# (Ideally after Step 7, 6, 11 are done)
# ============================================================
WORK_DIR="/path/to/project"
echo "Generating reads-only report..."
echo "Start time: $(date)"
python3 "${WORK_DIR}/scripts/19_reads_report.py" "$WORK_DIR"
echo "Finished: $(date)"
echo "Report: ${WORK_DIR}/processing/report/reads_report.html"

#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# ============================================================
# Step 17: Generate comprehensive HTML report
# Submit as: sbatch --dependency=afterok:<OVERVIEW>:<AMR_OVERVIEW> 17_generate_report.sh
# ============================================================
WORK_DIR="/path/to/project"
echo "Generating pipeline report..."
echo "Start time: $(date)"
python3 "${WORK_DIR}/scripts/17_generate_report_v2.py" "$WORK_DIR"
echo "Finished: $(date)"
echo "Report: ${WORK_DIR}/processing/report/pipeline_report.html"

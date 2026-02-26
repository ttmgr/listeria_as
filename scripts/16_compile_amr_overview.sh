#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# ============================================================
# Step 16: Compile AMRFinderPlus results into overview tables
# Submit as: sbatch --dependency=afterok:<AMRFINDER_JOB> 16_compile_amr_overview.sh
# ============================================================
WORK_DIR="/path/to/project"
echo "Compiling AMRFinderPlus overview..."
echo "Start time: $(date)"
python3 "${WORK_DIR}/scripts/16_compile_amr_overview.py" "$WORK_DIR"
echo "Finished: $(date)"

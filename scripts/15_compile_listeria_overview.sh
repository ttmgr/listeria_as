#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# -----------------------------------------------------------------------------
# Step 15: Build one integrated Listeria overview table and plots.
# Input: Listeria read summary, Listeria contig summary, and NanoStat outputs
# Output: processing/listeria/overview/listeria_overview.csv (+ figures)
# Run: sbatch --dependency=afterok:<LISTERIA_READS>:<LISTERIA_CONTIGS>:<NANOSTAT> scripts/15_compile_listeria_overview.sh
# -----------------------------------------------------------------------------
WORK_DIR="/path/to/project"
echo "Compiling comprehensive Listeria overview..."
echo "Start time: $(date)"
python3 "${WORK_DIR}/scripts/15_compile_listeria_overview.py" "$WORK_DIR"
echo "Finished: $(date)"

#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# ============================================================
# Plot Listeria results
# Submit as: sbatch 12_plot_listeria.sh
# ============================================================
WORK_DIR="/path/to/project"
INPUT="${WORK_DIR}/processing/listeria/listeria_summary.tsv"
OUTPUT="${WORK_DIR}/processing/listeria/plots"
echo "Plotting Listeria results..."
echo "Start time: $(date)"
python3 "${WORK_DIR}/scripts/plot_listeria.py" "$INPUT" "$OUTPUT"
echo "Finished: $(date)"

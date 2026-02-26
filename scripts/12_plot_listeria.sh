#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# -----------------------------------------------------------------------------
# Optional plot step: make quick Listeria read-count and base-count figures.
# Input: processing/listeria/listeria_summary.tsv
# Output: processing/listeria/plots/*.png and *.pdf
# Run: sbatch scripts/12_plot_listeria.sh
# -----------------------------------------------------------------------------
WORK_DIR="/path/to/project"
INPUT="${WORK_DIR}/processing/listeria/listeria_summary.tsv"
OUTPUT="${WORK_DIR}/processing/listeria/plots"
echo "Plotting Listeria results..."
echo "Start time: $(date)"
python3 "${WORK_DIR}/scripts/plot_listeria.py" "$INPUT" "$OUTPUT"
echo "Finished: $(date)"

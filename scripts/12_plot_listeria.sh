#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=120G
#SBATCH -t 00:30:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=plot_listeria
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Plot Listeria results
# Submit as: sbatch 12_plot_listeria.sh
# ============================================================

WORK_DIR="/lustre/groups/hpc/urban_lab/projects/tim"
INPUT="${WORK_DIR}/processing/listeria/listeria_summary.tsv"
OUTPUT="${WORK_DIR}/processing/listeria/plots"

echo "Plotting Listeria results..."
echo "Start time: $(date)"

python3 "${WORK_DIR}/scripts/plot_listeria.py" "$INPUT" "$OUTPUT"

echo "Finished: $(date)"

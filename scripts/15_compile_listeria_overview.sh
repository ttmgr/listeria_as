#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=120G
#SBATCH -t 1:00:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=listeria_overview
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 15: Compile comprehensive Listeria overview
# Submit as: sbatch --dependency=afterok:<LISTERIA_READS>:<LISTERIA_CONTIGS>:<NANOSTAT> 15_compile_listeria_overview.sh
# ============================================================

WORK_DIR="/lustre/groups/hpc/urban_lab/projects/tim"

echo "Compiling comprehensive Listeria overview..."
echo "Start time: $(date)"

python3 "${WORK_DIR}/scripts/15_compile_listeria_overview.py" "$WORK_DIR"

echo "Finished: $(date)"

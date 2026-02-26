#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=120G
#SBATCH -t 1:00:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=amr_overview
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 16: Compile AMRFinderPlus results into overview tables
# Submit as: sbatch --dependency=afterok:<AMRFINDER_JOB> 16_compile_amr_overview.sh
# ============================================================

WORK_DIR="/lustre/groups/hpc/urban_lab/projects/tim"

echo "Compiling AMRFinderPlus overview..."
echo "Start time: $(date)"

python3 "${WORK_DIR}/scripts/16_compile_amr_overview.py" "$WORK_DIR"

echo "Finished: $(date)"

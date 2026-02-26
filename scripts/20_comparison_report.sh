#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=120G
#SBATCH -t 1:00:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=comparison_report
#SBATCH -o /home/haicu/ttreska57/logs/%x_%j.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%j.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 20: Black Sample Comparison Report (AS vs N)
# Depends on: Steps 15, 16, 18 (all summary data computed)
# ============================================================

BASE_DIR="/lustre/groups/hpc/urban_lab/projects/tim"
SCRIPT_DIR="$(dirname "$0")"

echo "Generating comparison report..."
echo "Start time: $(date)"

python3 "${SCRIPT_DIR}/20_comparison_report.py" "$BASE_DIR"

echo "Finished: $(date)"

#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=120G
#SBATCH -t 0:30:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=reads_report
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 19: Reads-Only QC Report (Reads + Listeria + AMR)
# Submit as: sbatch 19_reads_report.sh
# (Ideally after Step 7, 6, 11 are done)
# ============================================================

WORK_DIR="/lustre/groups/hpc/urban_lab/projects/tim"

echo "Generating reads-only report..."
echo "Start time: $(date)"

python3 "${WORK_DIR}/scripts/19_reads_report.py" "$WORK_DIR"

echo "Finished: $(date)"
echo "Report: ${WORK_DIR}/processing/report/reads_report.html"

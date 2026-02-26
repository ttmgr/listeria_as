#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=120G
#SBATCH -t 0:45:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=generate_report_v2
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 17 (v2): Compile final HTML report
# Submit as: sbatch 17_generate_report_v2.sh
# ============================================================

WORK_DIR="/lustre/groups/hpc/urban_lab/projects/tim"

echo "Generating pipeline report (v2)..."
echo "Start time: $(date)"

# Run the v2 script
python3 "${WORK_DIR}/scripts/17_generate_report_v2.py" "$WORK_DIR"

echo "Finished: $(date)"
echo "Report: ${WORK_DIR}/processing/report/pipeline_report.html"

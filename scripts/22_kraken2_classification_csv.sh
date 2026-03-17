#!/bin/bash
# -----------------------------------------------------------------------------
# Step 22: Export Kraken2 read and contig classifications to CSV.
# Input: processing/kraken2/classified_*.txt and processing/kraken2_contigs/*/classified_*.txt
# Output: processing/kraken2_csv/{reads,flye,mdbg,myloasm}_classification.csv
# Run: sbatch --dependency=afterok:<KRAKEN_READS>:<KRAKEN_CONTIGS> scripts/22_kraken2_classification_csv.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

echo "Exporting Kraken2 classifications to CSV..."
echo "Start time: $(date)"

python3 "${SCRIPT_DIR}/22_kraken2_classification_csv.py" "$WORK_DIR"
EXIT_CODE=$?

if [ "$EXIT_CODE" -ne 0 ]; then
    ERRLOG="${SCRIPT_DIR}/error_22_kraken2_csv_$(date +%Y%m%d_%H%M%S).log"
    echo "ERROR: 22_kraken2_classification_csv.py exited with code ${EXIT_CODE}" > "$ERRLOG"
    echo "Timestamp: $(date)" >> "$ERRLOG"
    echo "Working dir: ${WORK_DIR}" >> "$ERRLOG"
    echo "SLURM Job ID: ${SLURM_JOB_ID:-local}" >> "$ERRLOG"
    echo "Error log written to: ${ERRLOG}"
    exit "$EXIT_CODE"
fi

echo "Finished: $(date)"

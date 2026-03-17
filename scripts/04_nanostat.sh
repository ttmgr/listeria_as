#!/bin/bash
# -----------------------------------------------------------------------------
# Step 4: Calculate per-sample read QC metrics (count, N50, quality, bases).
# Input: processing/nanofilt/filtered_<sample>.fastq
# Output: processing/nanostat/nanostat_<sample>.txt
# Run: sbatch --array=1-N scripts/04_nanostat.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

INPUT_DIR="${WORK_DIR}/processing/nanofilt"
OUTPUT_DIR="${WORK_DIR}/processing/nanostat"
mkdir -p "$OUTPUT_DIR"
# Get the filename for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(derive_basename "$BAM_FILE")
INPUT_FASTQ="${INPUT_DIR}/filtered_${BASENAME}.fastq"
if [ ! -f "$INPUT_FASTQ" ]; then
    echo "ERROR: Input file $INPUT_FASTQ does not exist. Skipping."
    exit 1
fi
echo "Processing: ${BASENAME}"
echo "Start time: $(date)"
NanoStat --fastq "$INPUT_FASTQ" \
         --name "${BASENAME}" \
         --outdir "$OUTPUT_DIR" \
         --tsv
# Rename output to include sample name
# NanoStat versions may produce different filenames: <name>NanoStats.txt or NanoStat-report.txt
RENAMED=0
for candidate in "${OUTPUT_DIR}/${BASENAME}NanoStats.txt" \
                 "${OUTPUT_DIR}/${BASENAME}NanoStat-report.txt" \
                 "${OUTPUT_DIR}/NanoStats.txt" \
                 "${OUTPUT_DIR}/NanoStat-report.txt" \
                 "${OUTPUT_DIR}/${BASENAME}"; do
    if [ -f "$candidate" ]; then
        mv "$candidate" "${OUTPUT_DIR}/nanostat_${BASENAME}.txt"
        echo "Renamed: $(basename "$candidate") -> nanostat_${BASENAME}.txt"
        RENAMED=1
        break
    fi
done
if [ "$RENAMED" -eq 0 ]; then
    echo "WARNING: No NanoStat output file found for ${BASENAME}. Checked:"
    echo "  ${BASENAME}NanoStats.txt, ${BASENAME}NanoStat-report.txt, NanoStats.txt, NanoStat-report.txt, ${BASENAME} (bare)"
    ls -la "${OUTPUT_DIR}/"*"${BASENAME}"* 2>/dev/null || echo "  No files matching ${BASENAME} in output dir."
fi
echo "Finished: $(date)"

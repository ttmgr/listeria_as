#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# ============================================================
# Step 6: Extract Listeria reads from Kraken2 output + stats
# Submit as: sbatch --array=1-66 --dependency=afterok:<KRAKEN_JOB_ID> 06_listeria_extract.sh
# ============================================================
INPUT_FASTQ_DIR="/path/to/project/processing/nanofilt"
KRAKEN_DIR="/path/to/project/processing/kraken2"
OUTPUT_DIR="/path/to/project/processing/listeria"
FILELIST="/path/to/project/filelist.txt"
mkdir -p "$OUTPUT_DIR"
# Get the filename for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(basename "$BAM_FILE" .bam)
INPUT_FASTQ="${INPUT_FASTQ_DIR}/filtered_${BASENAME}.fastq"
KRAKEN_OUTPUT="${KRAKEN_DIR}/classified_${BASENAME}.txt"
KRAKEN_REPORT="${KRAKEN_DIR}/report_${BASENAME}.txt"
if [ ! -f "$KRAKEN_OUTPUT" ]; then
    echo "ERROR: Kraken2 output $KRAKEN_OUTPUT does not exist. Skipping."
    exit 1
fi
echo "Processing: ${BASENAME}"
echo "Start time: $(date)"
# ---- Extract read IDs classified under Listeria (genus level) ----
# Kraken2 output format: C/U \t readID \t taxon \t length \t kmer_info
# Filter for lines containing "Listeria" in the taxon field (column 3)
grep -i "Listeria" "$KRAKEN_OUTPUT" | awk '{print $2}' > "${OUTPUT_DIR}/listeria_readids_${BASENAME}.txt"
LISTERIA_COUNT=$(wc -l < "${OUTPUT_DIR}/listeria_readids_${BASENAME}.txt")
echo "Listeria reads found: $LISTERIA_COUNT"
if [ "$LISTERIA_COUNT" -gt 0 ]; then
    # Extract Listeria reads from the FASTQ using seqtk
    seqtk subseq "$INPUT_FASTQ" "${OUTPUT_DIR}/listeria_readids_${BASENAME}.txt" \
        > "${OUTPUT_DIR}/listeria_${BASENAME}.fastq"
    # Run NanoStat on Listeria reads
    NanoStat --fastq "${OUTPUT_DIR}/listeria_${BASENAME}.fastq" \
             --name "listeria_${BASENAME}" \
             --outdir "$OUTPUT_DIR" \
             --tsv
    # Compute quick summary including median
    LISTERIA_BASES=$(awk 'NR%4==2 {sum += length($0)} END {print sum}' "${OUTPUT_DIR}/listeria_${BASENAME}.fastq")
    LISTERIA_MEAN_LEN=$(awk 'NR%4==2 {sum += length($0); n++} END {if(n>0) printf "%.1f", sum/n; else print 0}' "${OUTPUT_DIR}/listeria_${BASENAME}.fastq")
    LISTERIA_MEDIAN_LEN=$(awk 'NR%4==2 {lens[NR] = length($0); n++} END {
        asort(lens);
        if(n==0) {print 0}
        else if(n%2==1) {print lens[int(n/2)+1]}
        else {printf "%.1f", (lens[n/2] + lens[n/2+1]) / 2}
    }' "${OUTPUT_DIR}/listeria_${BASENAME}.fastq")
    echo "Listeria reads: $LISTERIA_COUNT"
    echo "Listeria bases: $LISTERIA_BASES"
    echo "Listeria mean read length: $LISTERIA_MEAN_LEN"
    echo "Listeria median read length: $LISTERIA_MEDIAN_LEN"
    # Atomic append to shared summary TSV (flock prevents interleaving)
    SUMMARY="${OUTPUT_DIR}/listeria_summary.tsv"
    (
        flock -x 200
        echo -e "${BASENAME}\t${LISTERIA_COUNT}\t${LISTERIA_BASES}\t${LISTERIA_MEAN_LEN}\t${LISTERIA_MEDIAN_LEN}" >> "$SUMMARY"
    ) 200>"${SUMMARY}.lock"
else
    echo "No Listeria reads found in ${BASENAME}"
    SUMMARY="${OUTPUT_DIR}/listeria_summary.tsv"
    (
        flock -x 200
        echo -e "${BASENAME}\t0\t0\t0\t0" >> "$SUMMARY"
    ) 200>"${SUMMARY}.lock"
fi
# ---- Also extract Listeria lines from Kraken2 report ----
grep -i "Listeria" "$KRAKEN_REPORT" > "${OUTPUT_DIR}/listeria_report_${BASENAME}.txt" 2>/dev/null || true
echo "Finished: $(date)"

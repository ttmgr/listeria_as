#!/bin/bash
# -----------------------------------------------------------------------------
# Step 14: Extract contigs classified as Listeria and summarize contig metrics.
# Input: Kraken2 contig classification files and assembly FASTA files
# Output: Listeria-only contigs plus processing/listeria/listeria_contigs_summary.tsv
# Run: sbatch --array=1-N --dependency=afterok:<KRAKEN_CONTIGS> scripts/14_listeria_contigs.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

KRAKEN_CONTIG_DIR="${WORK_DIR}/processing/kraken2_contigs"
FLYE_DIR="${WORK_DIR}/processing/racon"
MDBG_DIR="${WORK_DIR}/processing/mdbg"
MYLOASM_DIR="${WORK_DIR}/processing/myloasm"
OUTPUT_DIR="${WORK_DIR}/processing/listeria"
SUMMARY_FILE="${OUTPUT_DIR}/listeria_contigs_summary.tsv"
mkdir -p "${OUTPUT_DIR}/contigs_flye" "${OUTPUT_DIR}/contigs_mdbg" "${OUTPUT_DIR}/contigs_myloasm"
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(derive_basename "$BAM_FILE")
echo "Processing: ${BASENAME}"
echo "Start time: $(date)"
# Function to extract Listeria contigs and compute stats
extract_listeria_contigs() {
    local source=$1
    local classification=$2
    local asm_fasta=$3
    local out_fasta=$4
    echo "[$source] Extracting Listeria contigs..."
    if [ -s "$out_fasta" ]; then
        echo "  Already extracted. Skipping (found $out_fasta)."
        return
    fi
    if [ ! -s "$classification" ]; then
        echo "  No Kraken2 file for ${source}. Skipping."
        (
            flock -x 200
            echo -e "${BASENAME}\t${source}\t0\t0\t0\t0" >> "$SUMMARY_FILE"
        ) 200>"${SUMMARY_FILE}.lock"
        return
    fi
    # Extract contig IDs classified as Listeria
    local ID_FILE="${OUTPUT_DIR}/contigs_${source}/listeria_contigids_${BASENAME}.txt"
    grep -i "Listeria" "$classification" | awk '{print $2}' > "$ID_FILE"
    local CONTIG_COUNT=$(wc -l < "$ID_FILE")
    echo "  ${source}: ${CONTIG_COUNT} Listeria contigs"
    if [ "$CONTIG_COUNT" -gt 0 ] && [ -f "$asm_fasta" ]; then
        # Extract Listeria contigs using awk (FASTA format)
        # We use split() on the FASTA header to extract just the first word (the ID)
        awk 'NR==FNR{ids[$1]; next} /^>/{split($1, a, " "); name=substr(a[1], 2); p=(name in ids)} p' \
            "$ID_FILE" "$asm_fasta" > "$out_fasta"
        # Compute contig stats — accumulate sequence per contig (handles wrapped FASTA)
        local TOTAL_BASES=$(awk '!/^>/{sum += length($0)} END {print sum+0}' "$out_fasta")
        local MEDIAN_LEN=$(awk '
            /^>/ { if (seq != "") lens[++n] = length(seq); seq = ""; next }
            { seq = seq $0 }
            END {
                if (seq != "") lens[++n] = length(seq)
                if (n == 0) { print 0; exit }
                # Simple insertion sort for small arrays
                for (i = 2; i <= n; i++) {
                    key = lens[i]; j = i - 1
                    while (j >= 1 && lens[j] > key) { lens[j+1] = lens[j]; j-- }
                    lens[j+1] = key
                }
                if (n % 2 == 1) print lens[int(n/2)+1]
                else printf "%.0f\n", (lens[n/2] + lens[n/2+1]) / 2
            }
        ' "$out_fasta")
        local TOTAL_CONTIGS=$(grep -c "^>" "$out_fasta" || echo 0)
        echo "  ${source}: ${TOTAL_CONTIGS} contigs, ${TOTAL_BASES} bases, median ${MEDIAN_LEN} bp"
        (
            flock -x 200
            echo -e "${BASENAME}\t${source}\t${TOTAL_CONTIGS}\t${TOTAL_BASES}\t${MEDIAN_LEN}\t${CONTIG_COUNT}" >> "$SUMMARY_FILE"
        ) 200>"${SUMMARY_FILE}.lock"
    else
        (
            flock -x 200
            echo -e "${BASENAME}\t${source}\t0\t0\t0\t0" >> "$SUMMARY_FILE"
        ) 200>"${SUMMARY_FILE}.lock"
    fi
}
# ---- Flye ----
extract_listeria_contigs "flye" \
    "${KRAKEN_CONTIG_DIR}/flye/classified_flye_${BASENAME}.txt" \
    "${FLYE_DIR}/polished_${BASENAME}.fasta" \
    "${OUTPUT_DIR}/contigs_flye/listeria_flye_${BASENAME}.fasta"
# ---- metaMDBG ----
extract_listeria_contigs "mdbg" \
    "${KRAKEN_CONTIG_DIR}/mdbg/classified_mdbg_${BASENAME}.txt" \
    "${MDBG_DIR}/${BASENAME}/contigs.fasta" \
    "${OUTPUT_DIR}/contigs_mdbg/listeria_mdbg_${BASENAME}.fasta"
# ---- Myloasm ----
extract_listeria_contigs "myloasm" \
    "${KRAKEN_CONTIG_DIR}/myloasm/classified_myloasm_${BASENAME}.txt" \
    "${MYLOASM_DIR}/${BASENAME}/assembly_primary.fa" \
    "${OUTPUT_DIR}/contigs_myloasm/listeria_myloasm_${BASENAME}.fasta"
echo "Finished: $(date)"

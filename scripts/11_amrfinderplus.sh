#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# ============================================================
# Step 11: AMRFinderPlus --plus on reads, Flye contigs, mdbg contigs, myloasm contigs
# Submit as: sbatch --array=1-66 --dependency=afterok:<SEQKIT>:<FLYE>:<MDBG>:<MYLOASM> 11_amrfinderplus.sh
# ============================================================
READS_DIR="/path/to/project/processing/fasta"
FLYE_DIR="/path/to/project/processing/racon"
MDBG_DIR="/path/to/project/processing/mdbg"
MYLOASM_DIR="/path/to/project/processing/myloasm"
OUTPUT_DIR="/path/to/project/processing/amrfinder"
FILELIST="/path/to/project/filelist.txt"
THREADS=8
mkdir -p "${OUTPUT_DIR}/reads" "${OUTPUT_DIR}/flye" "${OUTPUT_DIR}/mdbg" "${OUTPUT_DIR}/myloasm"
# Get the filename for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(basename "$BAM_FILE" .bam)
echo "Processing: ${BASENAME}"
echo "Start time: $(date)"
# ---- 1. AMRFinderPlus on reads (FASTA) ----
READS_FASTA="${READS_DIR}/${BASENAME}.fasta"
if [ -f "$READS_FASTA" ]; then
    if [ -s "${OUTPUT_DIR}/reads/amrfinder_reads_${BASENAME}.tsv" ]; then
        echo "[Reads @ $(date)] Already processed. Skipping."
    else
        echo "[Reads @ $(date)] Running AMRFinderPlus --plus on reads..."
        amrfinder --plus \
            -n "$READS_FASTA" \
            --threads $THREADS \
            -o "${OUTPUT_DIR}/reads/amrfinder_reads_${BASENAME}.tsv"
        echo "  Done. $(wc -l < ${OUTPUT_DIR}/reads/amrfinder_reads_${BASENAME}.tsv) hits"
    fi
else
    echo "WARNING: Reads FASTA not found: $READS_FASTA"
fi
# ---- 2. AMRFinderPlus on Flye+Racon polished contigs ----
FLYE_CONTIGS="${FLYE_DIR}/polished_${BASENAME}.fasta"
if [ -f "$FLYE_CONTIGS" ]; then
    if [ -s "${OUTPUT_DIR}/flye/amrfinder_flye_${BASENAME}.tsv" ]; then
        echo "[Flye @ $(date)] Already processed. Skipping."
    else
        echo "[Flye @ $(date)] Running AMRFinderPlus --plus on Flye contigs..."
        amrfinder --plus \
            -n "$FLYE_CONTIGS" \
            --threads $THREADS \
            -o "${OUTPUT_DIR}/flye/amrfinder_flye_${BASENAME}.tsv"
        # Fallback if empty due to dummy FASTA
        if [ ! -s "${OUTPUT_DIR}/flye/amrfinder_flye_${BASENAME}.tsv" ]; then
            echo -e "Protein identifier\tContig id\tStart\tStop\tStrand\tElement symbol\tSequence name\tScope\tElement type\tElement subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence\tCov\tId\tAlign\tAcc\tName\tHMM\tDesc" > "${OUTPUT_DIR}/flye/amrfinder_flye_${BASENAME}.tsv"
        fi
        echo "  Done. $(wc -l < ${OUTPUT_DIR}/flye/amrfinder_flye_${BASENAME}.tsv) hits"
    fi
else
    echo "WARNING: Flye contigs not found: $FLYE_CONTIGS"
fi
# ---- 3. AMRFinderPlus on metaMDBG contigs ----
# metaMDBG outputs .fasta.gz — decompress if needed
MDBG_CONTIGS_GZ="${MDBG_DIR}/${BASENAME}/contigs.fasta.gz"
MDBG_CONTIGS="${MDBG_DIR}/${BASENAME}/contigs.fasta"
if [ -f "$MDBG_CONTIGS_GZ" ] && [ ! -f "$MDBG_CONTIGS" ]; then
    echo "Decompressing mdbg contigs..."
    gunzip -k "$MDBG_CONTIGS_GZ"
fi
if [ -f "$MDBG_CONTIGS" ]; then
    if [ -s "${OUTPUT_DIR}/mdbg/amrfinder_mdbg_${BASENAME}.tsv" ]; then
        echo "[MDBG @ $(date)] Already processed. Skipping."
    else
        echo "[MDBG @ $(date)] Running AMRFinderPlus --plus on mdbg contigs..."
        amrfinder --plus \
            -n "$MDBG_CONTIGS" \
            --threads $THREADS \
            -o "${OUTPUT_DIR}/mdbg/amrfinder_mdbg_${BASENAME}.tsv"
        # Fallback if empty due to dummy FASTA
        if [ ! -s "${OUTPUT_DIR}/mdbg/amrfinder_mdbg_${BASENAME}.tsv" ]; then
            echo -e "Protein identifier\tContig id\tStart\tStop\tStrand\tElement symbol\tSequence name\tScope\tElement type\tElement subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence\tCov\tId\tAlign\tAcc\tName\tHMM\tDesc" > "${OUTPUT_DIR}/mdbg/amrfinder_mdbg_${BASENAME}.tsv"
        fi
        echo "  Done. $(wc -l < ${OUTPUT_DIR}/mdbg/amrfinder_mdbg_${BASENAME}.tsv) hits"
    fi
else
    echo "WARNING: MDBG contigs not found for ${BASENAME}"
fi
# ---- 4. AMRFinderPlus on Myloasm contigs ----
MYLOASM_CONTIGS="${MYLOASM_DIR}/${BASENAME}/assembly_primary.fa"
if [ -f "$MYLOASM_CONTIGS" ]; then
    if [ -s "${OUTPUT_DIR}/myloasm/amrfinder_myloasm_${BASENAME}.tsv" ]; then
        echo "[Myloasm @ $(date)] Already processed. Skipping."
    else
        echo "[Myloasm @ $(date)] Running AMRFinderPlus --plus on Myloasm contigs..."
        amrfinder --plus \
            -n "$MYLOASM_CONTIGS" \
            --threads $THREADS \
            -o "${OUTPUT_DIR}/myloasm/amrfinder_myloasm_${BASENAME}.tsv"
        # Fallback if empty due to dummy FASTA
        if [ ! -s "${OUTPUT_DIR}/myloasm/amrfinder_myloasm_${BASENAME}.tsv" ]; then
            echo -e "Protein identifier\tContig id\tStart\tStop\tStrand\tElement symbol\tSequence name\tScope\tElement type\tElement subtype\tClass\tSubclass\tMethod\tTarget length\tReference sequence\tCov\tId\tAlign\tAcc\tName\tHMM\tDesc" > "${OUTPUT_DIR}/myloasm/amrfinder_myloasm_${BASENAME}.tsv"
        fi
        echo "  Done. $(wc -l < ${OUTPUT_DIR}/myloasm/amrfinder_myloasm_${BASENAME}.tsv) hits"
    fi
else
    echo "WARNING: Myloasm contigs not found for ${BASENAME}"
fi
echo "Finished all AMRFinderPlus analyses: $(date)"

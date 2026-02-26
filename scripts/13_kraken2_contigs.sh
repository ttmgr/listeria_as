#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=250G
#SBATCH -t 8:00:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=kraken2_contigs
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A_%a.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A_%a.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 13: Kraken2 on assembled contigs (Flye + mdbg + myloasm)
# Submit as: sbatch --array=1-66 --dependency=afterok:<FLYE>:<MDBG>:<MYLOASM> 13_kraken2_contigs.sh
# ============================================================


FLYE_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/racon"
MDBG_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/mdbg"
MYLOASM_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/myloasm"
KRAKEN_DB="/lustre/groups/hpc/urban_lab/datasets/ncbi/kraken2_core"
OUTPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/kraken2_contigs"
FILELIST="/lustre/groups/hpc/urban_lab/projects/tim/filelist.txt"
THREADS=20

mkdir -p "${OUTPUT_DIR}/flye" "${OUTPUT_DIR}/mdbg" "${OUTPUT_DIR}/myloasm"

BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(basename "$BAM_FILE" .bam)

echo "Processing: ${BASENAME}"
echo "Start time: $(date)"

# ---- 1. Kraken2 on Flye+Racon polished contigs ----
FLYE_CONTIGS="${FLYE_DIR}/polished_${BASENAME}.fasta"
if [ -s "$FLYE_CONTIGS" ]; then
    if [ -s "${OUTPUT_DIR}/flye/classified_flye_${BASENAME}.txt" ]; then
        echo "[Flye @ $(date)] Already classified. Skipping."
    else
        echo "[Flye @ $(date)] Running Kraken2 on Flye contigs..."
        kraken2 --db "$KRAKEN_DB" \
                --threads $THREADS \
                --use-names \
                --output "${OUTPUT_DIR}/flye/classified_flye_${BASENAME}.txt" \
                --report "${OUTPUT_DIR}/flye/report_flye_${BASENAME}.txt" \
                "$FLYE_CONTIGS"
        
        # Fallback if kraken2 outputs empty file due to empty/dummy FASTA
        if [ ! -s "${OUTPUT_DIR}/flye/classified_flye_${BASENAME}.txt" ]; then
            echo -e "U\tdummy\t0\t0\t0" > "${OUTPUT_DIR}/flye/classified_flye_${BASENAME}.txt"
            echo -e "100.00\t0\t0\tU\t0\tunclassified" > "${OUTPUT_DIR}/flye/report_flye_${BASENAME}.txt"
        fi
    
    if [ -f "${OUTPUT_DIR}/flye/classified_flye_${BASENAME}.txt" ]; then
        echo "  Done. $(wc -l < ${OUTPUT_DIR}/flye/classified_flye_${BASENAME}.txt) contigs classified"
    else
        echo "  Done. 0 contigs classified (Output empty/missing)."
    fi
    fi
else
    echo "WARNING: Flye contigs missing or empty: $FLYE_CONTIGS"
fi

# ---- 2. Kraken2 on metaMDBG contigs ----
# metaMDBG outputs .fasta.gz — decompress if needed
MDBG_CONTIGS_GZ="${MDBG_DIR}/${BASENAME}/contigs.fasta.gz"
MDBG_CONTIGS="${MDBG_DIR}/${BASENAME}/contigs.fasta"
if [ -s "$MDBG_CONTIGS_GZ" ] && [ ! -f "$MDBG_CONTIGS" ]; then
    echo "Decompressing mdbg contigs..."
    gunzip -k "$MDBG_CONTIGS_GZ"
fi

if [ -s "$MDBG_CONTIGS" ]; then
    if [ -s "${OUTPUT_DIR}/mdbg/classified_mdbg_${BASENAME}.txt" ]; then
        echo "[MDBG @ $(date)] Already classified. Skipping."
    else
        echo "[MDBG @ $(date)] Running Kraken2 on mdbg contigs..."
        kraken2 --db "$KRAKEN_DB" \
                --threads $THREADS \
                --use-names \
                --output "${OUTPUT_DIR}/mdbg/classified_mdbg_${BASENAME}.txt" \
                --report "${OUTPUT_DIR}/mdbg/report_mdbg_${BASENAME}.txt" \
                "$MDBG_CONTIGS"

        # Fallback if kraken2 outputs empty file due to empty/dummy FASTA
        if [ ! -s "${OUTPUT_DIR}/mdbg/classified_mdbg_${BASENAME}.txt" ]; then
            echo -e "U\tdummy\t0\t0\t0" > "${OUTPUT_DIR}/mdbg/classified_mdbg_${BASENAME}.txt"
            echo -e "100.00\t0\t0\tU\t0\tunclassified" > "${OUTPUT_DIR}/mdbg/report_mdbg_${BASENAME}.txt"
        fi

    if [ -f "${OUTPUT_DIR}/mdbg/classified_mdbg_${BASENAME}.txt" ]; then
        echo "  Done. $(wc -l < ${OUTPUT_DIR}/mdbg/classified_mdbg_${BASENAME}.txt) contigs classified"
    else
        echo "  Done. 0 contigs classified (Output empty/missing)."
    fi
    fi
else
    echo "WARNING: MDBG contigs missing or empty for ${BASENAME}"
fi

# ---- 3. Kraken2 on Myloasm contigs ----
MYLOASM_CONTIGS="${MYLOASM_DIR}/${BASENAME}/assembly_primary.fa"

if [ -s "$MYLOASM_CONTIGS" ]; then
    if [ -s "${OUTPUT_DIR}/myloasm/classified_myloasm_${BASENAME}.txt" ]; then
        echo "[Myloasm @ $(date)] Already classified. Skipping."
    else
        echo "[Myloasm @ $(date)] Running Kraken2 on Myloasm contigs..."
        kraken2 --db "$KRAKEN_DB" \
                --threads $THREADS \
                --use-names \
                --output "${OUTPUT_DIR}/myloasm/classified_myloasm_${BASENAME}.txt" \
                --report "${OUTPUT_DIR}/myloasm/report_myloasm_${BASENAME}.txt" \
                "$MYLOASM_CONTIGS"

        # Fallback if kraken2 outputs empty file due to empty/dummy FASTA
        if [ ! -s "${OUTPUT_DIR}/myloasm/classified_myloasm_${BASENAME}.txt" ]; then
            echo -e "U\tdummy\t0\t0\t0" > "${OUTPUT_DIR}/myloasm/classified_myloasm_${BASENAME}.txt"
            echo -e "100.00\t0\t0\tU\t0\tunclassified" > "${OUTPUT_DIR}/myloasm/report_myloasm_${BASENAME}.txt"
        fi

    if [ -f "${OUTPUT_DIR}/myloasm/classified_myloasm_${BASENAME}.txt" ]; then
        echo "  Done. $(wc -l < ${OUTPUT_DIR}/myloasm/classified_myloasm_${BASENAME}.txt) contigs classified"
    else
        echo "  Done. 0 contigs classified (Output empty/missing)."
    fi
    fi
else
    echo "WARNING: Myloasm contigs missing or empty for ${BASENAME}"
fi

echo "Finished: $(date)"

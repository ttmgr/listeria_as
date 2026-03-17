#!/bin/bash
# ============================================================
# Step 9b: Dorado Polish — assembly polishing with original reads
#
# This script aligns the original filtered reads back to a draft
# assembly using dorado aligner, then runs dorado polish to
# produce a higher-accuracy consensus assembly.
#
# It loops over all three assemblers (Flye, MetaMDBG, Myloasm)
# and polishes each draft if a contigs file exists.
#
# SLURM usage:
#   sbatch --array=1-66 scripts/09b_dorado_polish.sh
#
# Local usage:
#   export SLURM_ARRAY_TASK_ID=1
#   export SLURM_CPUS_PER_TASK=8
#   bash scripts/09b_dorado_polish.sh
# ============================================================

# --- Config ---
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

READS_DIR="${WORK_DIR}/processing/nanofilt"
THREADS="${SLURM_CPUS_PER_TASK:-8}"

# --- Dorado Models ---
export DORADO_MODELS_DIRECTORY="${DORADO_MODELS_DIR}"

# --- Resolve sample ID from filelist ---
BAM_PATH=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
ID=$(derive_basename "$BAM_PATH")

echo "=============================="
echo "  Dorado Polish: ${ID}"
echo "  Threads: ${THREADS}"
echo "=============================="

# Input reads (filtered FASTQ from step 3)
READS="${READS_DIR}/filtered_${ID}.fastq"

if [ ! -f "$READS" ]; then
    echo "ERROR: Filtered reads not found: $READS"
    exit 1
fi

# Output base directory for polished assemblies
POLISH_BASE="${WORK_DIR}/processing/assemblies"

# --- Loop over assemblers ---
for ASSEMBLER in flye mdbg myloasm; do

    echo ""
    echo "--- Assembler: ${ASSEMBLER} ---"

    # Locate the draft assembly (matching actual output paths from steps 8, 8b, 9)
    case "$ASSEMBLER" in
        flye)
            DRAFT="${WORK_DIR}/processing/racon/polished_${ID}.fasta"
            ;;
        mdbg)
            # metaMDBG outputs .fasta.gz — decompress if needed
            MDBG_GZ="${WORK_DIR}/processing/mdbg/${ID}/contigs.fasta.gz"
            DRAFT="${WORK_DIR}/processing/mdbg/${ID}/contigs.fasta"
            if [ -s "$MDBG_GZ" ] && [ ! -f "$DRAFT" ]; then
                echo "  Decompressing metaMDBG contigs..."
                gunzip -k "$MDBG_GZ"
            fi
            ;;
        myloasm)
            DRAFT="${WORK_DIR}/processing/myloasm/${ID}/assembly_primary.fa"
            ;;
    esac

    OUTDIR="${POLISH_BASE}/${ASSEMBLER}_polished/${ID}"
    ALIGNED="${OUTDIR}/aligned_reads.bam"
    POLISHED="${OUTDIR}/polished_assembly.fasta"

    # Skip if draft does not exist or is empty
    if [ ! -s "$DRAFT" ]; then
        echo "  No draft assembly found at ${DRAFT}, skipping."
        continue
    fi

    # Skip if already polished
    if [ -s "$POLISHED" ]; then
        echo "  Already polished: ${POLISHED}, skipping."
        continue
    fi

    mkdir -p "$OUTDIR"

    # Step 1: Align reads to draft with dorado aligner, sort and index
    echo "  Aligning reads to draft..."
    dorado aligner "$DRAFT" "$READS" \
        | samtools sort --threads "$THREADS" > "$ALIGNED"
    samtools index "$ALIGNED"

    # Step 2: Polish with dorado polish (--bacteria for bacterial genomes)
    echo "  Polishing assembly..."
    dorado polish "$ALIGNED" "$DRAFT" --bacteria > "$POLISHED"

    # Verify output
    if [ -s "$POLISHED" ]; then
        echo "  Done: ${POLISHED}"
    else
        echo "  WARNING: Polished output is empty."
    fi

done

echo ""
echo "Dorado polish complete for ${ID}."

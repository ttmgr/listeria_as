#!/bin/bash
# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim
# -----------------------------------------------------------------------------
# Step 18: Summarize assembly quality metrics with seqkit stats.
# Input: assembly FASTA files from Flye, metaMDBG, and Myloasm
# Output: processing/stats/assembly_stats_{flye,mdbg,myloasm}.tsv
# Run: sbatch scripts/18_assembly_stats.sh
# -----------------------------------------------------------------------------
PROJECT_DIR="/path/to/project"
MDBG_DIR="${PROJECT_DIR}/processing/mdbg"
FLYE_DIR="${PROJECT_DIR}/processing/flye"
MYLOASM_DIR="${PROJECT_DIR}/processing/myloasm"
STATS_DIR="${PROJECT_DIR}/processing/stats"
mkdir -p "$STATS_DIR"
echo "Calculating assembly statistics..."
echo "Start time: $(date)"
# 1. metaMDBG stats
# Find all contigs.fasta.gz files
echo "Processing metaMDBG..."
find "$MDBG_DIR" -name "contigs.fasta.gz" | sort > "${STATS_DIR}/mdbg_files.txt"
if [ -s "${STATS_DIR}/mdbg_files.txt" ]; then
    seqkit stats -a -T -j 4 $(cat "${STATS_DIR}/mdbg_files.txt") > "${STATS_DIR}/assembly_stats_mdbg.tsv"
    echo "Saved: ${STATS_DIR}/assembly_stats_mdbg.tsv"
else
    echo "No metaMDBG contigs found."
fi
# 2. Flye stats
# Find all assembly.fasta files
echo "Processing Flye..."
find "$FLYE_DIR" -name "assembly.fasta" | sort > "${STATS_DIR}/flye_files.txt"
if [ -s "${STATS_DIR}/flye_files.txt" ]; then
    seqkit stats -a -T -j 4 $(cat "${STATS_DIR}/flye_files.txt") > "${STATS_DIR}/assembly_stats_flye.tsv"
    echo "Saved: ${STATS_DIR}/assembly_stats_flye.tsv"
else
    echo "No Flye assemblies found."
fi
# 3. Myloasm stats
# Find all assembly.fasta files
echo "Processing Myloasm..."
find "$MYLOASM_DIR" -name "assembly_primary.fa" | sort > "${STATS_DIR}/myloasm_files.txt"
if [ -s "${STATS_DIR}/myloasm_files.txt" ]; then
    seqkit stats -a -T -j 4 $(cat "${STATS_DIR}/myloasm_files.txt") > "${STATS_DIR}/assembly_stats_myloasm.tsv"
    echo "Saved: ${STATS_DIR}/assembly_stats_myloasm.tsv"
else
    echo "No Myloasm assemblies found."
fi
echo "Finished: $(date)"

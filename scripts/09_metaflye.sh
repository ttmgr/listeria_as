#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=120G
#SBATCH -t 24:00:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=metaflye
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A_%a.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A_%a.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 9: metaFlye assembly + minimap2 alignment + racon polishing
# Submit as: sbatch --array=1-66 09_metaflye.sh
# ============================================================


INPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/nanofilt"
FLYE_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/flye"
MINIMAP_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/minimap2"
SAMTOOLS_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/samtools_bam"
RACON_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/racon"
FILELIST="/lustre/groups/hpc/urban_lab/projects/tim/filelist.txt"
THREADS=20

# Get the filename for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(basename "$BAM_FILE" .bam)
INPUT_FASTQ="${INPUT_DIR}/filtered_${BASENAME}.fastq"

if [ ! -f "$INPUT_FASTQ" ]; then
    echo "ERROR: Input file $INPUT_FASTQ does not exist. Skipping."
    exit 1
fi

mkdir -p "$FLYE_DIR" "$MINIMAP_DIR" "$SAMTOOLS_DIR" "$RACON_DIR"

FLYE_OUT="${FLYE_DIR}/assembly_${BASENAME}"
ASSEMBLY="${FLYE_OUT}/assembly.fasta"
SAM_FILE="${MINIMAP_DIR}/aligned_${BASENAME}.sam"
BAM_SORTED="${SAMTOOLS_DIR}/sorted_${BASENAME}.bam"
POLISHED="${RACON_DIR}/polished_${BASENAME}.fasta"

echo "Processing: ${BASENAME}"
echo "Start time: $(date)"

# Step 1: Flye meta-assembly
echo "[Step 1 @ $(date)] Running metaFlye..."
flye --meta --nano-hq "$INPUT_FASTQ" \
     --threads $THREADS \
     -o "$FLYE_OUT"

if [ ! -f "$ASSEMBLY" ]; then
    echo "Flye assembly aborted (likely low reads). Creating empty bypass files."
    echo ">no_assembly_low_reads_flye" > "$ASSEMBLY"
    echo ">no_assembly_low_reads_flye" > "$POLISHED"
else
    echo "Flye contigs: $(grep -c '>' ${ASSEMBLY} 2>/dev/null || echo 0)"

    # Step 2: minimap2 alignment
    echo "[Step 2 @ $(date)] Running minimap2..."
    minimap2 -ax map-ont -t $THREADS "$ASSEMBLY" "$INPUT_FASTQ" > "$SAM_FILE"
    
    # Step 3: SAM to sorted BAM
    echo "[Step 3 @ $(date)] Converting SAM → sorted BAM..."
    samtools view -b -@ $THREADS "$SAM_FILE" | samtools sort -@ $THREADS -o "$BAM_SORTED"
    samtools index "$BAM_SORTED"
    
    # Step 4: Racon polishing
    echo "[Step 4 @ $(date)] Running Racon polishing..."
    racon -t $THREADS "$INPUT_FASTQ" "$SAM_FILE" "$ASSEMBLY" > "$POLISHED"
    echo "Polished contigs: $(grep -c '>' ${POLISHED} 2>/dev/null || echo 0)"
    
    # Cleanup SAM (large file, BAM is kept)
    rm -f "$SAM_FILE"
fi

echo "Finished: $(date)"

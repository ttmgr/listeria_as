#!/bin/bash
#SBATCH -p cpu_p
#SBATCH -q cpu_normal
#SBATCH --mem=120G
#SBATCH -t 20:00:00
#SBATCH --nice=10000
#SBATCH --mail-user=timthilomaria.reska@helmholtz-munich.de
#SBATCH --mail-type=FAIL
#SBATCH -c 20
#SBATCH --job-name=metamdbg
#SBATCH -o /home/haicu/ttreska57/logs/%x_%A_%a.out
#SBATCH -e /home/haicu/ttreska57/logs/%x_%A_%a.err

# Activate conda environment
source ~/miniconda3/etc/profile.d/conda.sh
conda activate tim

# ============================================================
# Step 8: metaMDBG assembly
# Submit as: sbatch --array=1-66 08_metamdbg.sh
# ============================================================


INPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/nanofilt"
OUTPUT_DIR="/lustre/groups/hpc/urban_lab/projects/tim/processing/mdbg"
FILELIST="/lustre/groups/hpc/urban_lab/projects/tim/filelist.txt"

# Get the filename for this array task
BAM_FILE=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "$FILELIST")
BASENAME=$(basename "$BAM_FILE" .bam)
INPUT_FASTQ="${INPUT_DIR}/filtered_${BASENAME}.fastq"
SAMPLE_OUT="${OUTPUT_DIR}/${BASENAME}"

if [ ! -f "$INPUT_FASTQ" ]; then
    echo "ERROR: Input file $INPUT_FASTQ does not exist. Skipping."
    exit 1
fi

mkdir -p "$SAMPLE_OUT"

echo "Processing: ${BASENAME}"
echo "Start time: $(date)"

metaMDBG asm \
    --out-dir "$SAMPLE_OUT" \
    --in-ont "$INPUT_FASTQ" \
    --threads 20

# Fallback: if metaMDBG aborted due to lack of reads, create a dummy file
if [ ! -f "${SAMPLE_OUT}/contigs.fasta.gz" ]; then
    echo "Assembly failed or aborted (likely low reads). Creating empty bypass file."
    echo ">no_assembly_low_reads_mdbg" > "${SAMPLE_OUT}/contigs.fasta"
    gzip "${SAMPLE_OUT}/contigs.fasta"
fi

echo "Finished: $(date)"
echo "Contigs: $(grep -c '>' ${SAMPLE_OUT}/contigs.fasta 2>/dev/null || echo 0)"

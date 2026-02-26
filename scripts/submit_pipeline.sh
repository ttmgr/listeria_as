#!/bin/bash
# ============================================================
# Smart pipeline submission — skips completed steps
# Usage: cd /path/to/project && bash scripts/submit_pipeline.sh
# ============================================================
BAM_DIR="/path/to/BAM_files"
WORK_DIR="/path/to/project"
SCRIPT_DIR="${WORK_DIR}/scripts"
cd "$WORK_DIR"
# ---- Step 0: Create/update file list ----
echo "Creating file list from: $BAM_DIR"
# Only find original raw bam files (typically barcodeXX_XXX.bam), excluding any intermediate .sorted.bam files
ls "$BAM_DIR"/*.bam | grep -v '\.sorted\.bam$' | xargs -n1 basename | sort > filelist.txt
NUM_FILES=$(wc -l < filelist.txt)
echo "Found $NUM_FILES BAM files"
echo ""
# Create output directories
mkdir -p processing/{samtools,porechop,nanofilt,nanostat,kraken2,listeria,stats,mdbg,myloasm,flye,minimap2,samtools_bam,racon,fasta,amrfinder/reads,amrfinder/flye,amrfinder/mdbg,amrfinder/myloasm,kraken2_contigs/flye,kraken2_contigs/mdbg,kraken2_contigs/myloasm,listeria/contigs_flye,listeria/contigs_mdbg,listeria/contigs_myloasm,listeria/overview,report}
# ============================================================
# Helper: count how many samples are done for a given step
# Args: $1 = output_dir, $2 = filename_pattern (use BASENAME placeholder)
# Returns the number of MISSING output files
# ============================================================
check_step() {
    local out_dir="$1"
    local pattern="$2"
    local count=0
    while read -r bam; do
        local base=$(basename "$bam" .bam)
        local expected=$(echo "$pattern" | sed "s/BASENAME/${base}/g")
        if [ ! -s "${WORK_DIR}/${out_dir}/${expected}" ]; then
            ((count++))
        fi
    done < filelist.txt
    echo "$count"
}
# Build array of which task IDs still need processing for a step
get_pending_tasks() {
    local out_dir="$1"
    local pattern="$2"
    local pending=()
    local task_id=1
    while read -r bam; do
        local base=$(basename "$bam" .bam)
        local expected=$(echo "$pattern" | sed "s/BASENAME/${base}/g")
        if [ ! -s "${WORK_DIR}/${out_dir}/${expected}" ]; then
            pending+=("$task_id")
        fi
        ((task_id++))
    done < filelist.txt
    # Convert to SLURM array format: 1,3,5,7 etc.
    local IFS=','
    echo "${pending[*]}"
}
echo "========================================="
echo "  Checking pipeline status..."
echo "========================================="
echo ""
# Track job IDs for dependency chaining
JOB1="" JOB2="" JOB3="" JOB4="" JOB5="" JOB6="" JOB7="" JOB8="" JOB8B="" JOB9="" JOB10="" JOB11="" JOB13="" JOB14="" JOB15="" JOB16="" JOB18="" JOB17=""
# ---- Step 1: samtools BAM → FASTQ ----
MISSING_1=$(check_step "processing/samtools" "BASENAME.fastq")
if [ "$MISSING_1" -eq 0 ]; then
    echo "✓ Step 1 (samtools)     — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_1=$(get_pending_tasks "processing/samtools" "BASENAME.fastq")
    JOB1=$(sbatch --array=${PENDING_1} --parsable ${SCRIPT_DIR}/01_samtools_bam2fastq.sh)
    echo "→ Step 1 (samtools)     — submitting $MISSING_1 jobs: $JOB1"
fi
# ---- Step 2: Porechop ----
MISSING_2=$(check_step "processing/porechop" "trimmed_BASENAME.fastq")
if [ "$MISSING_2" -eq 0 ]; then
    echo "✓ Step 2 (porechop)     — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_2=$(get_pending_tasks "processing/porechop" "trimmed_BASENAME.fastq")
    DEP_2=""
    [ -n "$JOB1" ] && DEP_2="--dependency=afterok:${JOB1}"
    JOB2=$(sbatch --array=${PENDING_2} ${DEP_2} --parsable ${SCRIPT_DIR}/02_porechop.sh)
    echo "→ Step 2 (porechop)     — submitting $MISSING_2 jobs: $JOB2"
fi
# ---- Step 3: NanoFilt ----
MISSING_3=$(check_step "processing/nanofilt" "filtered_BASENAME.fastq")
if [ "$MISSING_3" -eq 0 ]; then
    echo "✓ Step 3 (nanofilt)     — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_3=$(get_pending_tasks "processing/nanofilt" "filtered_BASENAME.fastq")
    DEP_3=""
    [ -n "$JOB2" ] && DEP_3="--dependency=afterok:${JOB2}"
    JOB3=$(sbatch --array=${PENDING_3} ${DEP_3} --parsable ${SCRIPT_DIR}/03_nanofilt.sh)
    echo "→ Step 3 (nanofilt)     — submitting $MISSING_3 jobs: $JOB3"
fi
# ---- Determine dependency for post-nanofilt steps ----
# If nanofilt was submitted, depend on it. Otherwise no dependency needed.
POST_NANOFILT_DEP=""
[ -n "$JOB3" ] && POST_NANOFILT_DEP="--dependency=afterok:${JOB3}"
# ---- Step 3b: Read length extraction ----
if [ ! -s "processing/read_lengths_filtered_agg.tsv" ]; then
    JOB3B=$(sbatch --array=1-${NUM_FILES} ${POST_NANOFILT_DEP} --parsable ${SCRIPT_DIR}/03b_read_lengths.sh)
    echo "→ Step 3b (read len)    — submitting: $JOB3B"
else
    echo "✓ Step 3b (read len)    — DONE"
fi
# ---- Step 4: NanoStat ----
MISSING_4=$(check_step "processing/nanostat" "BASENAME")
if [ "$MISSING_4" -eq 0 ]; then
    echo "✓ Step 4 (nanostat)     — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_4=$(get_pending_tasks "processing/nanostat" "BASENAME")
    JOB4=$(sbatch --array=${PENDING_4} ${POST_NANOFILT_DEP} --parsable ${SCRIPT_DIR}/04_nanostat.sh)
    echo "→ Step 4 (nanostat)     — submitting $MISSING_4 jobs: $JOB4"
fi
# ---- Step 5: Kraken2 ----
MISSING_5=$(check_step "processing/kraken2" "report_BASENAME.txt")
if [ "$MISSING_5" -eq 0 ]; then
    echo "✓ Step 5 (kraken2)      — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_5=$(get_pending_tasks "processing/kraken2" "report_BASENAME.txt")
    JOB5=$(sbatch --array=${PENDING_5} ${POST_NANOFILT_DEP} --parsable ${SCRIPT_DIR}/05_kraken2.sh)
    echo "→ Step 5 (kraken2)      — submitting $MISSING_5 jobs: $JOB5"
fi
# ---- Step 6: Listeria extraction ----
MISSING_6=$(check_step "processing/listeria" "listeria_readids_BASENAME.txt")
if [ "$MISSING_6" -eq 0 ]; then
    echo "✓ Step 6 (listeria)     — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_6=$(get_pending_tasks "processing/listeria" "listeria_readids_BASENAME.txt")
    DEP_6=""
    [ -n "$JOB5" ] && DEP_6="--dependency=afterok:${JOB5}"
    JOB6=$(sbatch --array=${PENDING_6} ${DEP_6} --parsable ${SCRIPT_DIR}/06_listeria_extract.sh)
    echo "→ Step 6 (listeria)     — submitting $MISSING_6 jobs: $JOB6"
fi
# ---- Step 7: Compile stats ----
# Always re-run to pick up any new results
COMPILE_DEPS=""
[ -n "$JOB4" ] && COMPILE_DEPS="${JOB4}"
[ -n "$JOB6" ] && COMPILE_DEPS="${COMPILE_DEPS:+${COMPILE_DEPS}:}${JOB6}"
if [ -n "$COMPILE_DEPS" ]; then
    JOB7=$(sbatch --dependency=afterok:${COMPILE_DEPS} --parsable ${SCRIPT_DIR}/07_compile_stats.sh)
    echo "→ Step 7 (compile)      — submitted: $JOB7 (after $COMPILE_DEPS)"
else
    echo "✓ Step 7 (compile)      — nothing new to compile"
fi
# ---- Step 8: metaMDBG assembly ----
# Check for contigs.fasta.gz inside per-sample output dirs (metaMDBG outputs gzipped)
MISSING_8=$(check_step "processing/mdbg" "BASENAME/contigs.fasta.gz")
if [ "$MISSING_8" -eq 0 ]; then
    echo "✓ Step 8 (metaMDBG)     — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_8=$(get_pending_tasks "processing/mdbg" "BASENAME/contigs.fasta.gz")
    JOB8=$(sbatch --array=${PENDING_8} ${POST_NANOFILT_DEP} --parsable ${SCRIPT_DIR}/08_metamdbg.sh)
    echo "→ Step 8 (metaMDBG)     — submitting $MISSING_8 jobs: $JOB8"
fi
# ---- Step 8b: Myloasm assembly ----
MISSING_8B=$(check_step "processing/myloasm" "BASENAME/assembly_primary.fa")
if [ "$MISSING_8B" -eq 0 ]; then
    echo "✓ Step 8b(Myloasm)      — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_8B=$(get_pending_tasks "processing/myloasm" "BASENAME/assembly_primary.fa")
    JOB8B=$(sbatch --array=${PENDING_8B} ${POST_NANOFILT_DEP} --parsable ${SCRIPT_DIR}/08b_myloasm.sh)
    echo "→ Step 8b(Myloasm)      — submitting $MISSING_8B jobs: $JOB8B"
fi
# ---- Step 9: metaFlye + minimap2 + racon ----
# Check for final polished output
MISSING_9=$(check_step "processing/racon" "polished_BASENAME.fasta")
if [ "$MISSING_9" -eq 0 ]; then
    echo "✓ Step 9 (metaFlye)     — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_9=$(get_pending_tasks "processing/racon" "polished_BASENAME.fasta")
    JOB9=$(sbatch --array=${PENDING_9} ${POST_NANOFILT_DEP} --parsable ${SCRIPT_DIR}/09_metaflye.sh)
    echo "→ Step 9 (metaFlye)     — submitting $MISSING_9 jobs: $JOB9"
fi
# ---- Step 10: seqkit FASTQ → FASTA ----
MISSING_10=$(check_step "processing/fasta" "BASENAME.fasta")
if [ "$MISSING_10" -eq 0 ]; then
    echo "✓ Step 10 (seqkit)      — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_10=$(get_pending_tasks "processing/fasta" "BASENAME.fasta")
    JOB10=$(sbatch --array=${PENDING_10} ${POST_NANOFILT_DEP} --parsable ${SCRIPT_DIR}/10_seqkit_fq2fa.sh)
    echo "→ Step 10 (seqkit)      — submitting $MISSING_10 jobs: $JOB10"
fi
# ---- Step 11: AMRFinderPlus --plus (reads + flye + mdbg + myloasm) ----
# Depends on seqkit (reads), metaFlye (contigs), metaMDBG (contigs), Myloasm (contigs)
MISSING_11=$(check_step "processing/amrfinder/myloasm" "amrfinder_myloasm_BASENAME.tsv")
if [ "$MISSING_11" -eq 0 ]; then
    echo "✓ Step 11 (AMRFinder)   — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_11=$(get_pending_tasks "processing/amrfinder/myloasm" "amrfinder_myloasm_BASENAME.tsv")
    AMR_DEPS=""
    [ -n "$JOB10" ] && AMR_DEPS="${JOB10}"
    [ -n "$JOB8" ] && AMR_DEPS="${AMR_DEPS:+${AMR_DEPS}:}${JOB8}"
    [ -n "$JOB8B" ] && AMR_DEPS="${AMR_DEPS:+${AMR_DEPS}:}${JOB8B}"
    [ -n "$JOB9" ] && AMR_DEPS="${AMR_DEPS:+${AMR_DEPS}:}${JOB9}"
    DEP_11=""
    [ -n "$AMR_DEPS" ] && DEP_11="--dependency=afterok:${AMR_DEPS}"
    JOB11=$(sbatch --array=${PENDING_11} ${DEP_11} --parsable ${SCRIPT_DIR}/11_amrfinderplus.sh)
    echo "→ Step 11 (AMRFinder)   — submitting $MISSING_11 jobs: $JOB11"
fi
# ---- Step 13: Kraken2 on contigs (Flye + mdbg + myloasm) ----
MISSING_13=$(check_step "processing/kraken2_contigs/myloasm" "classified_myloasm_BASENAME.txt")
if [ "$MISSING_13" -eq 0 ]; then
    echo "✓ Step 13 (kraken2 ctg) — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_13=$(get_pending_tasks "processing/kraken2_contigs/myloasm" "classified_myloasm_BASENAME.txt")
    CTG_DEPS=""
    [ -n "$JOB8" ] && CTG_DEPS="${JOB8}"
    [ -n "$JOB8B" ] && CTG_DEPS="${CTG_DEPS:+${CTG_DEPS}:}${JOB8B}"
    [ -n "$JOB9" ] && CTG_DEPS="${CTG_DEPS:+${CTG_DEPS}:}${JOB9}"
    DEP_13=""
    [ -n "$CTG_DEPS" ] && DEP_13="--dependency=afterok:${CTG_DEPS}"
    JOB13=$(sbatch --array=${PENDING_13} ${DEP_13} --parsable ${SCRIPT_DIR}/13_kraken2_contigs.sh)
    echo "→ Step 13 (kraken2 ctg) — submitting $MISSING_13 jobs: $JOB13"
fi
# ---- Step 14: Extract Listeria contigs + stats ----
MISSING_14=$(check_step "processing/listeria/contigs_myloasm" "listeria_myloasm_BASENAME.fasta")
if [ "$MISSING_14" -eq 0 ]; then
    echo "✓ Step 14 (list. ctg)  — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_14=$(get_pending_tasks "processing/listeria/contigs_myloasm" "listeria_myloasm_BASENAME.fasta")
    DEP_14=""
    [ -n "$JOB13" ] && DEP_14="--dependency=afterok:${JOB13}"
    JOB14=$(sbatch --array=${PENDING_14} ${DEP_14} --parsable ${SCRIPT_DIR}/14_listeria_contigs.sh)
    echo "→ Step 14 (list. ctg)  — submitting $MISSING_14 jobs: $JOB14"
fi
# ---- Step 15: Comprehensive Listeria overview ----
# Depends on Listeria reads (6), Listeria contigs (14), NanoStat (4)
OVR_DEPS=""
[ -n "$JOB6" ] && OVR_DEPS="${JOB6}"
[ -n "$JOB14" ] && OVR_DEPS="${OVR_DEPS:+${OVR_DEPS}:}${JOB14}"
[ -n "$JOB4" ] && OVR_DEPS="${OVR_DEPS:+${OVR_DEPS}:}${JOB4}"
DEP_15=""
[ -n "$OVR_DEPS" ] && DEP_15="--dependency=afterok:${OVR_DEPS}"
JOB15=$(sbatch ${DEP_15} --parsable ${SCRIPT_DIR}/15_compile_listeria_overview.sh)
echo "→ Step 15 (overview)   — submitting: $JOB15"
# ---- Step 16: Compile AMRFinderPlus overview ----
DEP_16=""
[ -n "$JOB11" ] && DEP_16="--dependency=afterok:${JOB11}"
JOB16=$(sbatch ${DEP_16} --parsable ${SCRIPT_DIR}/16_compile_amr_overview.sh)
echo "→ Step 16 (AMR overview) — submitting: $JOB16"
echo ""
echo "========================================="
echo "  Pipeline dependency graph"
echo "========================================="
echo ""
# ---- Step 18: General Assembly Stats ----
STATS_MDBG="processing/stats/assembly_stats_mdbg.tsv"
STATS_FLYE="processing/stats/assembly_stats_flye.tsv"
STATS_MYLOASM="processing/stats/assembly_stats_myloasm.tsv"
if [ -f "$STATS_MDBG" ] && [ -f "$STATS_FLYE" ] && [ -f "$STATS_MYLOASM" ]; then
    echo "✓ Step 18 (Assembly Stats) — DONE"
else
    DEP_18=""
    [ -n "$JOB8" ] && DEP_18="${JOB8}"
    [ -n "$JOB8B" ] && DEP_18="${DEP_18:+${DEP_18}:}${JOB8B}"
    [ -n "$JOB9" ] && DEP_18="${DEP_18:+${DEP_18}:}${JOB9}"
    [ -n "$DEP_18" ] && DEP_18="--dependency=afterok:${DEP_18}"
    JOB18=$(sbatch ${DEP_18} --parsable ${SCRIPT_DIR}/18_assembly_stats.sh)
    echo "→ Step 18 (Assembly Stats) — submitting: $JOB18"
fi
# ---- Step 17: Generate HTML Report ----
DEP_17=""
# Depends on everything that feeds into the report
[ -n "$JOB15" ] && DEP_17="${JOB15}"
[ -n "$JOB16" ] && DEP_17="${DEP_17:+${DEP_17}:}${JOB16}"
[ -n "$JOB18" ] && DEP_17="${DEP_17:+${DEP_17}:}${JOB18}"
[ -n "$DEP_17" ] && DEP_17="--dependency=afterok:${DEP_17}"
JOB17=$(sbatch ${DEP_17} --parsable ${SCRIPT_DIR}/17_generate_report.sh)
echo "→ Step 17 (HTML Report)  — submitting: $JOB17"
# ---- Step 20: Black Sample Comparison Report ----
DEP_20=""
[ -n "$JOB17" ] && DEP_20="--dependency=afterok:${JOB17}"
JOB20=$(sbatch ${DEP_20} --parsable ${SCRIPT_DIR}/20_comparison_report.sh)
echo "→ Step 20 (Comparison)  — submitting: $JOB20"
echo ""
echo "========================================="
echo "  Pipeline dependency graph"
echo "========================================="
echo ""
echo "  1 samtools → 2 porechop → 3 nanofilt ─┬─ 4 nanostat ──┬─ 7 compile"
echo "                                         ├─ 5 kraken2 → 6 listeria ──────────┐"
echo "                                         ├─ 8 metaMDBG ──┬─ 13 kraken2_ctg ──┤"
echo "                                         ├─ 8b myloasm ──┤                   │"
echo "                                         ├─ 9 metaFlye ──┘      │            │"
echo "                                         │              14 list. contigs ─────┤"
echo "                                         └─ 10 seqkit ─┬─ 11 AMRFinder+ ─────┤"
echo "                                                     ├─ 15 list. overview ──┴─ 17 REPORT"
echo "                                                     ├─ 16 AMR overview ──────┘"
echo "                                                     └─ 18 assembly stats ────┘"
echo ""
echo "Monitor: squeue -u \$USER"

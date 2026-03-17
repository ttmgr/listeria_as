#!/bin/bash
# -----------------------------------------------------------------------------
# Main launcher for the Listeria adaptive-sampling workflow.
# What it does: submits all pipeline steps with dependencies and skips completed outputs.
# Before running: place BAM files into listeria_1/, listeria_2/, etc.
# Usage: cd <project_root> && bash scripts/submit_pipeline.sh [partition]
# -----------------------------------------------------------------------------
SCRIPT_DIR="$(cd "$(dirname "$0")" && pwd)"
export SCRIPT_DIR

# Capture user overrides BEFORE sourcing pipeline.conf (conda activate can reset env)
_SKIP_FILELIST="${SKIP_FILELIST:-0}"
_BARCODE_FILTER="${BARCODE_FILTER:-}"
_RUN_DORADO="${RUN_DORADO:-0}"

source "${SCRIPT_DIR}/pipeline.conf"

FILELIST_PATH="${WORK_DIR}/filelist.txt"
BARCODE_FILTER="$_BARCODE_FILTER"
RUN_DORADO="$_RUN_DORADO"
SKIP_FILELIST="$_SKIP_FILELIST"
LOG_DIR="/home/haicu/ttreska57/logs"

# Default SBATCH settings
SBATCH_ARGS=(
    --partition=cpu_p
    --qos=cpu_normal
    --nice=10000
    --mail-user=timthilomaria.reska@helmholtz-munich.de
    --mail-type=FAIL
    -o "${LOG_DIR}/%x_%A_%a.out"
    -e "${LOG_DIR}/%x_%A_%a.err"
)

# Standard resources: 120G / 12 cores / 24h
SBATCH_STD=(--mem=120G -c 12 -t 24:00:00)
# Heavy resources for Kraken2: 250G / 20 cores
SBATCH_HEAVY=(--mem=250G -c 20 -t 24:00:00)

cd "$WORK_DIR"
mkdir -p "$LOG_DIR"

# ---- Step 0b: Build sample_metadata.csv from master sample sheet ----
if [ -f "${SAMPLE_SHEET_DIR}/sample_sheet_master.csv" ]; then
    echo "Building sample_metadata.csv..."
    python3 "${SCRIPT_DIR}/00b_build_sample_metadata.py" "$WORK_DIR"
else
    echo "WARNING: No sample_sheet_master.csv found in ${SAMPLE_SHEET_DIR}"
    echo "  Cohort reports will not have metadata enrichment."
fi
echo ""

# ---- Step 0: Build file list from round directories ----
if [ "$SKIP_FILELIST" = "1" ] && [ -s "$FILELIST_PATH" ]; then
    echo "Using existing filelist (SKIP_FILELIST=1): $(wc -l < "$FILELIST_PATH") entries"
else
    echo "Building filelist from rounds: ${ROUNDS}"
    > "$FILELIST_PATH"
    for ROUND in $ROUNDS; do
        ROUND_DIR="${BAM_BASE_DIR}/listeria_${ROUND}"
        if [ ! -d "$ROUND_DIR" ]; then
            echo "WARNING: Round directory not found: $ROUND_DIR — skipping"
            continue
        fi
        find "$ROUND_DIR" -maxdepth 1 -type f -name '*.bam' ! -name '*.sorted.bam' | \
            sed "s#^${BAM_BASE_DIR}/##" | sort >> "$FILELIST_PATH"
    done

    if [ -n "$BARCODE_FILTER" ]; then
        echo "Filtering file list for: ${BARCODE_FILTER}"
        grep -i "$BARCODE_FILTER" "$FILELIST_PATH" > "${FILELIST_PATH}.tmp" || true
        mv "${FILELIST_PATH}.tmp" "$FILELIST_PATH"
    fi
fi

NUM_FILES=$(wc -l < "$FILELIST_PATH")
if [ "$NUM_FILES" -eq 0 ]; then
    echo "ERROR: No BAM files found for the current selection."
    exit 1
fi
echo "Found $NUM_FILES BAM files"
echo ""

# Create output directories
mkdir -p processing/{samtools,porechop,nanofilt,nanostat,kraken2,kraken2_csv,listeria,stats,mdbg,myloasm,flye,minimap2,samtools_bam,racon,fasta,dorado/flye,dorado/mdbg,dorado/myloasm,dorado/done,amrfinder/reads,amrfinder/flye,amrfinder/mdbg,amrfinder/myloasm,kraken2_contigs/flye,kraken2_contigs/mdbg,kraken2_contigs/myloasm,listeria/contigs_flye,listeria/contigs_mdbg,listeria/contigs_myloasm,listeria/overview,report}

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
        local base=$(derive_basename "$bam")
        local expected=$(echo "$pattern" | sed "s/BASENAME/${base}/g")
        if [ ! -s "${WORK_DIR}/${out_dir}/${expected}" ]; then
            ((count++))
        fi
    done < "$FILELIST_PATH"
    echo "$count"
}

# Build array of which task IDs still need processing for a step
get_pending_tasks() {
    local out_dir="$1"
    local pattern="$2"
    local pending=()
    local task_id=1
    while read -r bam; do
        local base=$(derive_basename "$bam")
        local expected=$(echo "$pattern" | sed "s/BASENAME/${base}/g")
        if [ ! -s "${WORK_DIR}/${out_dir}/${expected}" ]; then
            pending+=("$task_id")
        fi
        ((task_id++))
    done < "$FILELIST_PATH"
    # Convert to SLURM array format: 1,3,5,7 etc.
    local IFS=','
    echo "${pending[*]}"
}

normalize_nanostat_outputs() {
    local nano_dir="${WORK_DIR}/processing/nanostat"
    [ -d "$nano_dir" ] || return 0

    while read -r bam; do
        local base=$(derive_basename "$bam")
        local expected="${nano_dir}/nanostat_${base}.txt"
        [ -s "$expected" ] && continue

        for candidate in \
            "${nano_dir}/${base}NanoStats.txt" \
            "${nano_dir}/${base}NanoStats" \
            "${nano_dir}/${base}_NanoStats.txt" \
            "${nano_dir}/${base}_NanoStats" \
            "${nano_dir}/${base}.txt" \
            "${nano_dir}/${base}.tsv" \
            "${nano_dir}/${base}" \
            "${nano_dir}/nanostat_${base}"
        do
            if [ -s "$candidate" ]; then
                mv "$candidate" "$expected"
                break
            fi
        done
    done < "$FILELIST_PATH"
}

normalize_nanostat_outputs

echo "========================================="
echo "  Checking pipeline status..."
echo "========================================="
echo ""

# Track job IDs for dependency chaining
JOB1="" JOB2="" JOB3="" JOB3B="" JOB4="" JOB5="" JOB6="" JOB7="" JOB8="" JOB8B="" JOB9="" JOB9B="" JOB10="" JOB11="" JOB13="" JOB14="" JOB15="" JOB16="" JOB18="" JOB17=""

# ---- Step 1: samtools BAM → FASTQ ----
MISSING_1=$(check_step "processing/samtools" "BASENAME.fastq")
if [ "$MISSING_1" -eq 0 ]; then
    echo "✓ Step 1 (samtools)     — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_1=$(get_pending_tasks "processing/samtools" "BASENAME.fastq")
    JOB1=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=samtools --array=${PENDING_1} --parsable ${SCRIPT_DIR}/01_samtools_bam2fastq.sh)
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
    JOB2=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=porechop --array=${PENDING_2} ${DEP_2} --parsable ${SCRIPT_DIR}/02_porechop.sh)
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
    JOB3=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=nanofilt --array=${PENDING_3} ${DEP_3} --parsable ${SCRIPT_DIR}/03_nanofilt.sh)
    echo "→ Step 3 (nanofilt)     — submitting $MISSING_3 jobs: $JOB3"
fi

# ---- Determine dependency for post-nanofilt steps ----
POST_NANOFILT_DEP=""
[ -n "$JOB3" ] && POST_NANOFILT_DEP="--dependency=afterok:${JOB3}"

# ---- Step 3b: Read length extraction ----
if [ ! -s "processing/read_lengths_filtered_agg.tsv" ]; then
    JOB3B=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=readlen --array=1-${NUM_FILES} ${POST_NANOFILT_DEP} --parsable ${SCRIPT_DIR}/03b_read_lengths.sh)
    echo "→ Step 3b (read len)    — submitting: $JOB3B"
else
    echo "✓ Step 3b (read len)    — DONE"
fi

# ---- Step 4: NanoStat ----
MISSING_4=$(check_step "processing/nanostat" "nanostat_BASENAME.txt")
if [ "$MISSING_4" -eq 0 ]; then
    echo "✓ Step 4 (nanostat)     — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_4=$(get_pending_tasks "processing/nanostat" "nanostat_BASENAME.txt")
    JOB4=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=nanostat --array=${PENDING_4} ${POST_NANOFILT_DEP} --parsable ${SCRIPT_DIR}/04_nanostat.sh)
    echo "→ Step 4 (nanostat)     — submitting $MISSING_4 jobs: $JOB4"
fi

# ---- Step 5: Kraken2 ----
MISSING_5=$(check_step "processing/kraken2" "report_BASENAME.txt")
if [ "$MISSING_5" -eq 0 ]; then
    echo "✓ Step 5 (kraken2)      — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_5=$(get_pending_tasks "processing/kraken2" "report_BASENAME.txt")
    JOB5=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_HEAVY[@]}" --job-name=kraken2 --array=${PENDING_5} ${POST_NANOFILT_DEP} --parsable ${SCRIPT_DIR}/05_kraken2.sh)
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
    JOB6=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=listeria --array=${PENDING_6} ${DEP_6} --parsable ${SCRIPT_DIR}/06_listeria_extract.sh)
    echo "→ Step 6 (listeria)     — submitting $MISSING_6 jobs: $JOB6"
fi

# ---- Step 7: Compile stats ----
# Rerun if dependencies were resubmitted OR if the output CSV is missing/header-only
COMPILE_DEPS=""
[ -n "$JOB4" ] && COMPILE_DEPS="${JOB4}"
[ -n "$JOB6" ] && COMPILE_DEPS="${COMPILE_DEPS:+${COMPILE_DEPS}:}${JOB6}"
STATS_CSV="${WORK_DIR}/processing/stats/read_metrics_summary.csv"
STATS_LINES=0
[ -f "$STATS_CSV" ] && STATS_LINES=$(wc -l < "$STATS_CSV")
if [ -n "$COMPILE_DEPS" ]; then
    JOB7=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=compile --dependency=afterok:${COMPILE_DEPS} --parsable ${SCRIPT_DIR}/07_compile_stats.sh)
    echo "→ Step 7 (compile)      — submitted: $JOB7 (after $COMPILE_DEPS)"
elif [ "$STATS_LINES" -le 1 ]; then
    JOB7=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=compile --parsable ${SCRIPT_DIR}/07_compile_stats.sh)
    echo "→ Step 7 (compile)      — submitted: $JOB7 (output CSV empty/missing, recompiling)"
else
    echo "✓ Step 7 (compile)      — DONE ($((STATS_LINES - 1)) samples in CSV)"
fi

# ---- Step 8: metaMDBG assembly ----
MISSING_8=$(check_step "processing/mdbg" "BASENAME/contigs.fasta.gz")
if [ "$MISSING_8" -eq 0 ]; then
    echo "✓ Step 8 (metaMDBG)     — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_8=$(get_pending_tasks "processing/mdbg" "BASENAME/contigs.fasta.gz")
    JOB8=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=metamdbg --array=${PENDING_8} ${POST_NANOFILT_DEP} --parsable ${SCRIPT_DIR}/08_metamdbg.sh)
    echo "→ Step 8 (metaMDBG)     — submitting $MISSING_8 jobs: $JOB8"
fi

# ---- Step 8b: Myloasm assembly ----
MISSING_8B=$(check_step "processing/myloasm" "BASENAME/assembly_primary.fa")
if [ "$MISSING_8B" -eq 0 ]; then
    echo "✓ Step 8b(Myloasm)      — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_8B=$(get_pending_tasks "processing/myloasm" "BASENAME/assembly_primary.fa")
    JOB8B=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=myloasm --array=${PENDING_8B} ${POST_NANOFILT_DEP} --parsable ${SCRIPT_DIR}/08b_myloasm.sh)
    echo "→ Step 8b(Myloasm)      — submitting $MISSING_8B jobs: $JOB8B"
fi

# ---- Step 9: metaFlye + minimap2 + racon ----
MISSING_9=$(check_step "processing/racon" "polished_BASENAME.fasta")
if [ "$MISSING_9" -eq 0 ]; then
    echo "✓ Step 9 (metaFlye)     — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_9=$(get_pending_tasks "processing/racon" "polished_BASENAME.fasta")
    JOB9=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=metaflye --array=${PENDING_9} ${POST_NANOFILT_DEP} --parsable ${SCRIPT_DIR}/09_metaflye.sh)
    echo "→ Step 9 (metaFlye)     — submitting $MISSING_9 jobs: $JOB9"
fi

# ---- Step 9b: Dorado Polish ----
if [ "$RUN_DORADO" != "1" ]; then
    echo "○ Step 9b(Dorado)       — skipped (RUN_DORADO=0)"
else
    MISSING_9B=$(check_step "processing/dorado/done" "BASENAME.done")
    if [ "$MISSING_9B" -eq 0 ]; then
        echo "✓ Step 9b(Dorado)       — DONE ($NUM_FILES/$NUM_FILES)"
    else
        PENDING_9B=$(get_pending_tasks "processing/dorado/done" "BASENAME.done")
        DORADO_DEPS=""
        [ -n "$JOB8" ] && DORADO_DEPS="${JOB8}"
        [ -n "$JOB8B" ] && DORADO_DEPS="${DORADO_DEPS:+${DORADO_DEPS}:}${JOB8B}"
        [ -n "$JOB9" ] && DORADO_DEPS="${DORADO_DEPS:+${DORADO_DEPS}:}${JOB9}"
        [ -n "$JOB3" ] && DORADO_DEPS="${DORADO_DEPS:+${DORADO_DEPS}:}${JOB3}"
        DEP_9B=""
        [ -n "$DORADO_DEPS" ] && DEP_9B="--dependency=afterok:${DORADO_DEPS}"
        JOB9B=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=dorado --array=${PENDING_9B} ${DEP_9B} --parsable ${SCRIPT_DIR}/09b_dorado_polish.sh)
        echo "→ Step 9b(Dorado)       — submitting $MISSING_9B jobs: $JOB9B"
    fi
fi

# ---- Step 10: seqkit FASTQ → FASTA ----
MISSING_10=$(check_step "processing/fasta" "BASENAME.fasta")
if [ "$MISSING_10" -eq 0 ]; then
    echo "✓ Step 10 (seqkit)      — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_10=$(get_pending_tasks "processing/fasta" "BASENAME.fasta")
    JOB10=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=seqkit --array=${PENDING_10} ${POST_NANOFILT_DEP} --parsable ${SCRIPT_DIR}/10_seqkit_fq2fa.sh)
    echo "→ Step 10 (seqkit)      — submitting $MISSING_10 jobs: $JOB10"
fi

# ---- Step 11: AMRFinderPlus --plus (reads + flye + mdbg + myloasm) ----
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
    JOB11=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=amrfinder --array=${PENDING_11} ${DEP_11} --parsable ${SCRIPT_DIR}/11_amrfinderplus.sh)
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
    JOB13=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_HEAVY[@]}" --job-name=kraken2ctg --array=${PENDING_13} ${DEP_13} --parsable ${SCRIPT_DIR}/13_kraken2_contigs.sh)
    echo "→ Step 13 (kraken2 ctg) — submitting $MISSING_13 jobs: $JOB13"
fi

# ---- Step 22: Kraken2 Classification CSV Export ----
DEP_22=""
[ -n "$JOB5" ] && DEP_22="${JOB5}"
[ -n "$JOB13" ] && DEP_22="${DEP_22:+${DEP_22}:}${JOB13}"
[ -n "$DEP_22" ] && DEP_22="--dependency=afterok:${DEP_22}"
JOB22=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=kraken2csv ${DEP_22} --parsable ${SCRIPT_DIR}/22_kraken2_classification_csv.sh)
echo "→ Step 22 (kraken2 csv) — submitting: $JOB22"

# ---- Step 14: Extract Listeria contigs + stats ----
MISSING_14=$(check_step "processing/listeria/contigs_myloasm" "listeria_myloasm_BASENAME.fasta")
if [ "$MISSING_14" -eq 0 ]; then
    echo "✓ Step 14 (list. ctg)  — DONE ($NUM_FILES/$NUM_FILES)"
else
    PENDING_14=$(get_pending_tasks "processing/listeria/contigs_myloasm" "listeria_myloasm_BASENAME.fasta")
    DEP_14=""
    [ -n "$JOB13" ] && DEP_14="--dependency=afterok:${JOB13}"
    JOB14=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=listctg --array=${PENDING_14} ${DEP_14} --parsable ${SCRIPT_DIR}/14_listeria_contigs.sh)
    echo "→ Step 14 (list. ctg)  — submitting $MISSING_14 jobs: $JOB14"
fi

# ---- Step 15: Comprehensive Listeria overview ----
OVR_DEPS=""
[ -n "$JOB6" ] && OVR_DEPS="${JOB6}"
[ -n "$JOB14" ] && OVR_DEPS="${OVR_DEPS:+${OVR_DEPS}:}${JOB14}"
[ -n "$JOB4" ] && OVR_DEPS="${OVR_DEPS:+${OVR_DEPS}:}${JOB4}"
DEP_15=""
[ -n "$OVR_DEPS" ] && DEP_15="--dependency=afterok:${OVR_DEPS}"
JOB15=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=listovr ${DEP_15} --parsable ${SCRIPT_DIR}/15_compile_listeria_overview.sh)
echo "→ Step 15 (overview)   — submitting: $JOB15"

# ---- Step 16: Compile AMRFinderPlus overview ----
DEP_16=""
[ -n "$JOB11" ] && DEP_16="--dependency=afterok:${JOB11}"
JOB16=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=amrovr ${DEP_16} --parsable ${SCRIPT_DIR}/16_compile_amr_overview.sh)
echo "→ Step 16 (AMR overview) — submitting: $JOB16"
echo ""

echo "========================================="
echo "  Assembly & Report Steps"
echo "========================================="
echo ""

# ---- Step 18: General Assembly Stats ----
DEP_18=""
[ -n "$JOB8" ] && DEP_18="${JOB8}"
[ -n "$JOB8B" ] && DEP_18="${DEP_18:+${DEP_18}:}${JOB8B}"
[ -n "$JOB9" ] && DEP_18="${DEP_18:+${DEP_18}:}${JOB9}"
[ -n "$DEP_18" ] && DEP_18="--dependency=afterok:${DEP_18}"
JOB18=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=asmstats ${DEP_18} --parsable ${SCRIPT_DIR}/18_assembly_stats.sh)
echo "→ Step 18 (Assembly Stats) — submitting: $JOB18"

# ---- Step 17: Generate HTML Report ----
DEP_17=""
[ -n "$JOB7" ] && DEP_17="${JOB7}"
[ -n "$JOB15" ] && DEP_17="${DEP_17:+${DEP_17}:}${JOB15}"
[ -n "$JOB16" ] && DEP_17="${DEP_17:+${DEP_17}:}${JOB16}"
[ -n "$JOB18" ] && DEP_17="${DEP_17:+${DEP_17}:}${JOB18}"
[ -n "$DEP_17" ] && DEP_17="--dependency=afterok:${DEP_17}"
JOB17=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=report ${DEP_17} --parsable ${SCRIPT_DIR}/17_generate_report.sh)
echo "→ Step 17 (HTML Report)  — submitting: $JOB17"

# ---- Step 20: Comparison Report ----
DEP_20=""
[ -n "$JOB17" ] && DEP_20="--dependency=afterok:${JOB17}"
JOB20=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=comparison ${DEP_20} --parsable ${SCRIPT_DIR}/20_comparison_report.sh)
echo "→ Step 20 (Comparison)  — submitting: $JOB20"

# ---- Step 23: Cohort Reports & Post-Processing ----
DEP_23=""
[ -n "$JOB17" ] && DEP_23="${JOB17}"
[ -n "$JOB22" ] && DEP_23="${DEP_23:+${DEP_23}:}${JOB22}"
[ -n "$JOB20" ] && DEP_23="${DEP_23:+${DEP_23}:}${JOB20}"
[ -n "$DEP_23" ] && DEP_23="--dependency=afterok:${DEP_23}"
JOB23=$(sbatch "${SBATCH_ARGS[@]}" "${SBATCH_STD[@]}" --job-name=cohort_rpt ${DEP_23} --parsable ${SCRIPT_DIR}/23_cohort_reports.sh)
echo "→ Step 23 (Cohort Rpts) — submitting: $JOB23"

echo ""
echo "========================================="
echo "  Pipeline dependency graph"
echo "========================================="
echo ""
echo "  1 samtools → 2 porechop → 3 nanofilt ─┬─ 4 nanostat ──┬─ 7 compile"
echo "                                         ├─ 5 kraken2 → 6 listeria ──────────┐"
echo "                                         │              22 kraken2 csv ──────┤"
echo "                                         ├─ 8 metaMDBG ──┬─ 13 kraken2_ctg ──┤"
echo "                                         ├─ 8b myloasm ──┤                   │"
echo "                                         ├─ 9 metaFlye ──┘      │            │"
echo "                                         │   └─ 9b Dorado       │            │"
echo "                                         │              14 list. contigs ─────┤"
echo "                                         └─ 10 seqkit ─┬─ 11 AMRFinder+ ─────┤"
echo "                                                     ├─ 15 list. overview ──┤"
echo "                                                     ├─ 16 AMR overview ────┴─ 17 REPORT ─┬─ 20 Comparison"
echo "                                                     └─ 18 assembly stats ───┘            └─ 23 Cohort Reports"
echo ""
echo "Rounds processed: ${ROUNDS}"
echo "Total BAM files: ${NUM_FILES}"
echo "Monitor: squeue -u \$USER"

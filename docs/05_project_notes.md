# Project Notes & Lessons Learned

## Project Overview

Listeria Adaptive Sampling Pipeline — Oxford Nanopore sequencing for detecting
*Listeria monocytogenes* from food safety samples. Compares Adaptive Sampling (AS)
vs Normal (N) sequencing across multiple extraction methods (Sponge, Cotton, Zymo)
and two cohorts (Black, Blue).

**Cluster**: HPC at Helmholtz Munich (`hpc-build01`)
**Project dir**: `/lustre/groups/hpc/urban_lab/projects/tim/foodsafety`
**Conda env**: `tim`
**Rounds**: Round 1 (`listeria_1/`) and Round 2 (`listeria_2/`) — different biological samples, NOT re-sequencing

---

## Round 1 / Round 2 Architecture

Rounds are separated by **directory**, not by filename:
- `listeria_1/barcode03_AS.bam` (round 1)
- `listeria_2/barcode03_AS.bam` (round 2) — **different biological sample**, same barcode number

`derive_basename()` in `pipeline.conf` produces round-prefixed names:
- `listeria_1/barcode03_AS.bam` -> `r1_barcode03_AS`
- `listeria_2/barcode06_N.bam` -> `r2_barcode06_N`

All downstream files use this `r{N}_barcode{NN}_{AS|N}` format.

---

## Mistakes & Fixes (Do Not Repeat)

### 1. Leading whitespace in `filelist.txt`

**What happened**: When creating `filelist.txt` manually using a heredoc (`cat > filelist.txt << 'EOF'`), leading spaces from terminal indentation were included in each line. This caused `derive_basename` to produce names like `r  listeria_1_barcode03_AS` (with spaces), creating broken output files.

**Fix**: Use `printf` or `sed` to ensure no leading whitespace:
```bash
# GOOD — guaranteed no leading whitespace
printf '%s\n' \
  "listeria_1/barcode03_AS.bam" \
  "listeria_1/barcode03_N.bam" > filelist.txt

# GOOD — strip whitespace after writing
sed -i 's/^[[:space:]]*//' filelist.txt

# VERIFY — cat -A shows raw content, no spaces before text
cat -A filelist.txt
```

**How to detect**: Output files have spaces in names: `r  listeria_1_barcode01_AS.fastq`

### 2. SLURM jobs from previous runs not cancelled

**What happened**: After fixing the filelist and resubmitting, old broken jobs from the previous submission were still running (status `DependencyNeverSatisfied` or `Dependency`). These old jobs created new empty/broken output files that interfered with the new run.

**Fix**: Always verify ALL old jobs are cancelled before resubmitting:
```bash
scancel -u $USER
squeue --me              # must show ZERO jobs
# only then:
SKIP_FILELIST=1 bash scripts/submit_pipeline.sh
```

### 3. `SKIP_FILELIST` env var lost after `conda activate`

**What happened**: `SKIP_FILELIST=1 bash scripts/submit_pipeline.sh` didn't work because `source pipeline.conf` runs `conda activate`, which can reset environment variables in some shell configurations.

**Fix**: `submit_pipeline.sh` now captures env vars *before* sourcing `pipeline.conf`:
```bash
_SKIP_FILELIST="${SKIP_FILELIST:-0}"
source "${SCRIPT_DIR}/pipeline.conf"
SKIP_FILELIST="$_SKIP_FILELIST"
```

### 4. NanoStat file naming mismatch

**What happened**: `04_nanostat.sh` produced files like `{BASENAME}NanoStats.txt` but `07_compile_stats.sh` only searched for `nanostat_*.txt`. The `mv` to rename was silently suppressed with `2>/dev/null || true`, so no error was visible but stats were missing.

**Fix**:
- `04_nanostat.sh`: Now tries multiple candidate filenames for the rename
- `07_compile_stats.sh`: Falls back to `*NanoStats.txt` pattern if `nanostat_*.txt` finds nothing
- `submit_pipeline.sh`: Has `normalize_nanostat_outputs()` that fixes naming before checking status

### 5. Hardcoded `ylim=15000` in report plots

**What happened**: The Listeria reads bar chart in `17_generate_report_v2.py` had `ylim=15000`, which cut off bars when actual values were higher.

**Fix**: Compute dynamically: `ylim = int(max_val * 1.2)`

### 6. Guppy references in report (never used Guppy)

**What happened**: The report had a methods box mentioning Guppy basecalling and tried to run `guppy_basecaller --version`. Pipeline uses Dorado, not Guppy.

**Fix**: Removed the entire software versions section and methods box from the report.

### 7. No cohort distinction in cluster report

**What happened**: The cluster report (`17_generate_report_v2.py`) pooled all samples together with no way to tell Black from Blue.

**Fix**:
- Report now loads `sample_metadata.csv` and adds Cohort/Group columns to all tables
- Generates per-cohort per-round reports: `pipeline_report_black_r1.html`, etc.

### 8. "Local" scripts couldn't run on cluster

**What happened**: `build_local_black_report.py`, `22_local_plots.py`, and the split scripts were designed for local use only. Running the full analysis required downloading data and running Spyder locally.

**Fix**:
- Created `23_cohort_reports.sh` wrapper that runs all post-processing scripts
- Added as Step 23 in `submit_pipeline.sh` with proper SLURM dependencies
- Made `build_local_black_report.py` robust (try/except for missing files)

### 9. Metadata column name mismatch between cluster and local

**What happened**: `00b_build_sample_metadata.py` output columns like `label_colour`, `kit_used` but report scripts expected `cohort`, `kit`, `group`.

**Fix**: Updated `00b_build_sample_metadata.py` to output the canonical column names that all scripts expect: `cohort`, `group`, `swab_type`, `kit`, `basename`, `round`, `barcode_label`.

### 10. Regex patterns didn't handle round-prefixed sample names

**What happened**: Patterns like `r"barcode(\d+)_(AS|N)"` didn't match `r1_barcode03_AS`.

**Fix**: Updated all regex patterns across scripts to optionally match the round prefix:
```python
r"(?:r\d+_)?barcode(\d+)_(AS|N)"
```

Files updated: `build_local_black_report.py`, `15_compile_listeria_overview.py`, `17_generate_report_v2.py`, `22_local_plots.py`

---

## Quick Reference: Running the Pipeline

### Full run (all samples, all rounds)
```bash
cd /lustre/groups/hpc/urban_lab/projects/tim/foodsafety
bash scripts/submit_pipeline.sh
```

### Specific barcodes only
```bash
# Write filelist manually (NO leading spaces!)
printf '%s\n' \
  "listeria_1/barcode03_AS.bam" \
  "listeria_1/barcode03_N.bam" \
  "listeria_2/barcode06_AS.bam" \
  "listeria_2/barcode06_N.bam" > filelist.txt

SKIP_FILELIST=1 bash scripts/submit_pipeline.sh
```

### Monitor
```bash
watch squeue --me
```

### Check logs on failure
```bash
# SLURM logs
ls ~/logs/<stepname>_<jobid>_*.{out,err}

# e.g. for porechop job 34614956, task 1:
cat ~/logs/porechop_34614956_1.err
cat ~/logs/porechop_34614956_1.out
```

### Clean restart
```bash
scancel -u $USER          # cancel everything
squeue --me               # verify zero jobs
rm -rf processing/        # nuke outputs
bash scripts/submit_pipeline.sh   # fresh start
```

---

## Pipeline Output Structure

```
processing/
  report/
    pipeline_report.html              # combined (all rounds, all cohorts)
    pipeline_report_r1.html           # round 1 only
    pipeline_report_r2.html           # round 2 only
    pipeline_report_black_r1.html     # Black cohort, round 1
    pipeline_report_blue_r1.html      # Blue cohort, round 1
    pipeline_report_black_r2.html     # Black cohort, round 2
    pipeline_report_blue_r2.html      # Blue cohort, round 2
    qc_metrics.csv
    listeria_reads_summary.csv
    assembly_stats_{flye,mdbg,myloasm}.csv
    amr_genes_{reads,contigs}.csv

local_black_report_r1/                # detailed Black cohort report (round 1)
local_black_report_r2/                # detailed Black cohort report (round 2)
local_blue_report_r1/                 # detailed Blue cohort report (round 1)
local_blue_report_r2/                 # detailed Blue cohort report (round 2)

plots_black/                          # publication figures for Black
plots_blue/                           # publication figures for Blue
```

---

## Adding New Rounds (3, 4, 5, ...)

1. Place BAMs in `listeria_3/`, `listeria_4/`, etc.
2. Add rows to `samplesheets/sample_sheet_master.csv` with `Round=3`, etc.
3. Update `pipeline.conf`: `ROUNDS="1 2 3 4 5"`
4. Run `bash scripts/submit_pipeline.sh` — it auto-discovers everything

No code changes needed. Reports will auto-generate per-round splits.

---

## Cohorts & Barcode Mapping

Barcode-to-cohort mapping is defined in `samplesheets/sample_sheet_master.csv`.
The same barcode number in different rounds = **different biological sample**.

`00b_build_sample_metadata.py` reads the master sheet and produces `sample_metadata.csv`
with all the columns needed by report scripts.

Cohorts in the data: **Black**, **Blue**, Control, Gillian.
Only Black and Blue are used in cohort reports. Gillian and Control appear only in the combined report.

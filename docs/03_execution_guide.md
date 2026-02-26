# Execution Guide: Running the Pipeline

Before running any script, you must edit the hardcoded placeholder paths. 

**Required Edits:**
Check `scripts/submit_pipeline.sh` and each `*.sh` script you intend to run. Replace placeholders like:
- `/path/to/project`
- `/path/to/logs`
- `/path/to/BAM_files`

*(Note: Most failed runs come from path mismatches!)*

---

## Approach 1: High-Performance Cluster (SLURM) — *Recommended*

Assembly and Kraken2 steps require substantial RAM and CPU processing. Full datasets are therefore heavily recommended to be executed via a SLURM orchestrator.

### SLURM Submission Flags (Used by `submit_pipeline.sh`)
- `--array=...`: run one script instance per sample/task index (e.g. `--array=1-66`).
- `--dependency=afterok:<jobid>`: start only after another job finishes successfully.
- `--parsable`: print job ID only (easy for dependency chaining).

### Running the Orchestrator
The `submit_pipeline.sh` file manages all dependencies and generates a `filelist.txt` automatically from your inputs.

```bash
bash scripts/submit_pipeline.sh
```
This script will submit jobs 1-20 in order. Depending on HPC loading, this can take a few hours to a day.

**To check on your jobs:**
```bash
squeue -u "$USER"
```

---

## Approach 2: Local / Manual Execution

For single-sample tests or local runs on small data sets (or servers lacking SLURM scheduler overhead), you can execute the bash scripts manually. 

Because the scripts internally rely on SLURM array task IDs, you must set these variables manually in your terminal before running the script:

```bash
export SLURM_ARRAY_TASK_ID=1
export SLURM_CPUS_PER_TASK=4
bash scripts/01_samtools_bam2fastq.sh
```

---

## Output Locations to Check First
After a successful run, check the following endpoints to confirm completion:
- `processing/report/pipeline_report.html` (The final comprehensive report)
- `processing/listeria/overview/listeria_overview.csv`
- `processing/amrfinder/overview/amr_reads_overview.csv`
- `processing/stats/read_metrics_summary.csv`

---

## Troubleshooting & Common Issues

- **`No file found for task ID`**
  Check `filelist.txt` exists and that `SLURM_ARRAY_TASK_ID` does not exceed the line count.
- **Missing Kraken2 output**
  Verify your `KRAKEN2_DB` variable path inside the script matches a valid built DB.
- **Empty assembly outputs**
  Low-read samples may trigger the generation of fallback "dummy" files (e.g. `empty.fasta`) to keep downstream scripts from crashing. This is intentional.
- **Report missing sections**
  Run our input checker script to identify which specific component failed:
  ```bash
  bash scripts/check_report_inputs.sh /path/to/project
  ```

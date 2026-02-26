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

## Approach 2: Standalone Workstation (Multi-Core Linux/macOS)

If you don't have a SLURM cluster but *do* have a powerful workstation (e.g. 32+ cores, 128GB+ RAM), you can run the pipeline locally using a simple `bash` loop.

Because the scripts internally rely on the `$SLURM_ARRAY_TASK_ID` variable to know which sample to process, we can simply trick the scripts into thinking they are running on a cluster by exporting that variable in a loop.

### Creating a Local Run Script
Save the following code as `run_local.sh` in the root of the repository:

```bash
#!/bin/bash
# run_local.sh - Executes the pipeline on a standalone workstation

# 1. Edit your paths in the scripts as usual
# 2. Make sure filelist.txt is generated (has one BAM file per line)
TOTAL_SAMPLES=$(wc -l < "/path/to/project/filelist.txt")

# Set how many cores each job should use
export SLURM_CPUS_PER_TASK=8

echo "Starting local execution for $TOTAL_SAMPLES samples..."

# Loop through every sample sequentially
for i in $(seq 1 $TOTAL_SAMPLES); do
    echo "Processing Task ID: $i"
    export SLURM_ARRAY_TASK_ID=$i
    
    # Run the steps you need (uncomment as necessary)
    bash scripts/01_samtools_bam2fastq.sh
    # bash scripts/02_porechop.sh
    # bash scripts/03_nanofilt.sh
    bash scripts/05_kraken2.sh
    bash scripts/06_listeria_extract.sh
    # bash scripts/09_metaflye.sh
done

echo "Sample processing complete. Now running aggregate steps..."

# Run the summary scripts (these only need to run once at the end)
bash scripts/15_compile_listeria_overview.sh
bash scripts/16_compile_qc_metrics.sh
bash scripts/17_generate_report.sh

echo "Pipeline Finished!"
```

Run this script in your terminal:
```bash
bash run_local.sh
```

*Note: This runs samples sequentially (one after another). If you want to run samples in parallel on a massive server, look into using `GNU parallel` instead of a `for` loop.*

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

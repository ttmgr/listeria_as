# Listeria Adaptive Sampling Pipeline

This repository contains a full analysis workflow for Oxford Nanopore sequencing data, focused on *Listeria* detection and characterization from mixed samples.

The scripts are written for HPC batch execution (SLURM), but you can also run them manually for small tests.

## Why this matters for food safety

*Listeria monocytogenes* is a foodborne pathogen that can survive in food-processing environments and can cause severe disease, especially in pregnant people, older adults, and immunocompromised patients. Rapidly detecting *Listeria* signal in sequencing data helps food safety teams act faster during contamination checks and outbreak investigations.

## What adaptive sampling is (Nanopore, plain language)

Adaptive sampling is a real-time enrichment approach during Nanopore sequencing. As DNA starts passing through a pore, the instrument basecalls the first part of the read, compares it to a target reference, and then decides:

- keep sequencing the read if it looks like target DNA
- eject the read early if it looks off-target

In practice, this can increase sequencing yield for targets of interest (here, *Listeria*) without physically enriching DNA in the wet lab first.

## Who this README is for

This guide is written for wet-lab users who may not work with command-line pipelines every day. It is focused on:

- what to edit before first run
- which command to run
- what files should appear when things work

## Quick start (recommended path)

1. Clone the repo and enter it.
2. Edit path variables in scripts (details below).
3. Confirm tools are available in your conda environment.
4. Run the orchestrator script:

```bash
bash scripts/submit_pipeline.sh
```

This orchestrator submits and links the full workflow (steps 1 to 20), and skips outputs that already exist.

## Before first run: required edits

Most scripts use placeholder paths like:

- `/path/to/project`
- `/path/to/BAM_files`

At minimum, update these in:

- `scripts/submit_pipeline.sh`
- all step scripts you plan to run (for example `scripts/01_samtools_bam2fastq.sh`, `scripts/05_kraken2.sh`, `scripts/17_generate_report.sh`)

The scripts also assume:

- a conda environment named `tim`
- `conda.sh` at `~/miniconda3/etc/profile.d/conda.sh`

If your environment name or conda location differs, update the first lines of each shell script.

## Required input files

- Raw Nanopore `.bam` files in one directory (set as `BAM_DIR` in `submit_pipeline.sh`)
- File naming that keeps sample/barcode identity (the scripts derive sample IDs from BAM names)

When `submit_pipeline.sh` runs, it creates `filelist.txt` automatically from `BAM_DIR`.

## Pipeline map (what each script does)

1. `01_samtools_bam2fastq.sh`: BAM to FASTQ
2. `02_porechop.sh`: adapter trimming
3. `03_nanofilt.sh`: length filtering (`<100 bp` removed)
4. `03b_read_lengths.sh`: read-length distributions
5. `04_nanostat.sh`: QC metrics per sample
6. `05_kraken2.sh`: taxonomic classification on reads
7. `06_listeria_extract.sh`: extract *Listeria* reads and per-sample summaries
8. `07_compile_stats.sh`: compile read and *Listeria* summary tables
9. `08_metamdbg.sh`, `08b_myloasm.sh`, `09_metaflye.sh`: assembly workflows
10. `10_seqkit_fq2fa.sh`: FASTQ to FASTA
11. `11_amrfinderplus.sh`: AMR calls on reads and contigs
12. `13_kraken2_contigs.sh`: taxonomy on assembled contigs
13. `14_listeria_contigs.sh`: extract *Listeria* contigs and contig stats
14. `15_compile_listeria_overview.sh`: integrated *Listeria* overview
15. `16_compile_amr_overview.sh`: AMR overview tables
16. `18_assembly_stats.sh`: assembly summary stats
17. `17_generate_report.sh`: HTML report
18. `20_comparison_report.sh`: AS vs N comparison report

Optional downstream scripts:

- `19_reads_report.sh`: read-focused quick report
- `21_statistical_analysis.py`: statistical testing
- `22_local_plots.py`: local publication-style figure generation
- `20_export_tables_to_xlsx.py`, `21_kraken2_to_spreadsheets.py`: spreadsheet exports

## Output locations to check first

After a successful run, start here:

- `processing/report/pipeline_report.html`
- `processing/listeria/overview/listeria_overview.csv`
- `processing/amrfinder/overview/amr_reads_overview.csv`
- `processing/amrfinder/overview/amr_contigs_overview.csv`
- `processing/stats/read_metrics_summary.csv`

## Running on SLURM vs local machine

### SLURM (recommended for full dataset)

Use:

```bash
bash scripts/submit_pipeline.sh
```

Monitor jobs:

```bash
squeue -u "$USER"
```

### Local/manual (small tests only)

Some scripts expect SLURM variables. For single-sample tests, set them manually:

```bash
export SLURM_ARRAY_TASK_ID=1
export SLURM_CPUS_PER_TASK=4
bash scripts/01_samtools_bam2fastq.sh
```

Assembly and Kraken2 steps can require substantial RAM/CPU, so full runs are better on HPC.

## Common issues and fast checks

- `No file found for task ID`: check `filelist.txt` and `SLURM_ARRAY_TASK_ID`
- missing Kraken2 output: verify `KRAKEN2_DB` path
- empty assembly outputs: low-read samples may trigger fallback dummy files by design
- report missing sections: run input checker

```bash
bash scripts/check_report_inputs.sh /path/to/project
```

## Minimal software list

Core tools used across scripts:

- samtools
- porechop
- NanoFilt
- NanoStat
- kraken2
- seqtk
- seqkit
- flye
- metaMDBG
- myloasm
- minimap2
- racon
- amrfinder
- Python 3 with pandas/numpy/matplotlib/scipy

## Final note

If you only change one thing before running: make sure every placeholder path is replaced with real paths for your system. Most failed runs come from path mismatches.

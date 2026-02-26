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

## Install everything (copy/paste)

This pipeline is designed for Linux/macOS with conda-compatible package management. For Windows, use WSL2.

### 1) Set up Bioconda channels (one time)

```bash
conda config --add channels bioconda
conda config --add channels conda-forge
conda config --set channel_priority strict
```

### 2) Create one environment with all required tools

```bash
conda create -n listeria_as \
  python=3.10 \
  samtools porechop nanofilt nanostat kraken2 seqtk seqkit \
  flye metamdbg myloasm minimap2 racon ncbi-amrfinderplus \
  pandas numpy scipy matplotlib \
  -c conda-forge -c bioconda --strict-channel-priority
```

Then activate:

```bash
conda activate listeria_as
```

If you want to keep using your existing `tim` environment, install the same package set into that environment instead of creating `listeria_as`.

### 3) Quick check that tools are available

```bash
samtools --version
kraken2 --version
flye --version
amrfinder --version
```

## Where to get each tool (official package/docs links)

- Bioconda setup: <https://bioconda.github.io/>
- samtools: <https://bioconda.github.io/recipes/samtools/README.html>
- porechop: <https://bioconda.github.io/recipes/porechop/README.html>
- NanoFilt: <https://bioconda.github.io/recipes/nanofilt/README.html>
- NanoStat: <https://bioconda.github.io/recipes/nanostat/README.html>
- kraken2: <https://bioconda.github.io/recipes/kraken2/README.html>
- seqtk: <https://bioconda.github.io/recipes/seqtk/README.html>
- seqkit: <https://bioconda.github.io/recipes/seqkit/README.html>
- Flye: <https://bioconda.github.io/recipes/flye/README.html>
- metaMDBG: <https://bioconda.github.io/recipes/metamdbg/README.html>
- Myloasm: <https://bioconda.github.io/recipes/myloasm/README.html>
- minimap2: <https://bioconda.github.io/recipes/minimap2/README.html>
- racon: <https://bioconda.github.io/recipes/racon/README.html>
- AMRFinderPlus: <https://bioconda.github.io/recipes/ncbi-amrfinderplus/README.html>

## Databases: where to get them and how to configure

### Kraken2 database (required for steps 5 and 13)

Current scripts use a `KRAKEN2_DB` (or `KRAKEN_DB`) variable in:

- `scripts/05_kraken2.sh`
- `scripts/13_kraken2_contigs.sh`

Set those variables to your Kraken2 database directory.

Official options:

1. Build standard Kraken2 DB from NCBI/RefSeq (official manual):

```bash
kraken2-build --standard --threads 24 --db /path/to/kraken2_standard
```

2. Use official prebuilt Kraken2/Bracken indexes (AWS index listed by Kraken2 wiki):

- <https://benlangmead.github.io/aws-indexes/k2>

Reference docs:

- Kraken2 manual: <https://github.com/DerrickWood/kraken2/wiki/Manual>

### AMRFinderPlus database (required for step 11)

Install/update the AMRFinderPlus database with:

```bash
amrfinder --update
```

Alternative (custom location):

```bash
amrfinder_update -d /path/to/amrfinder_db
```

Then use:

```bash
amrfinder -d /path/to/amrfinder_db/latest ...
```

Reference docs:

- AMRFinder homepage: <https://www.ncbi.nlm.nih.gov/pathogens/antimicrobial-resistance/AMRFinder/>
- AMRFinder wiki: <https://github.com/ncbi/amr/wiki>

## Flag reference (all core pipeline command flags)

This section explains the command-line flags used in the core workflow scripts.

### SLURM submission flags (used by `submit_pipeline.sh`)

- `--array=...`: run one script instance per sample/task index.
- `--dependency=afterok:<jobid>`: start only after another job finishes successfully.
- `--parsable`: print job ID only (easy for dependency chaining).

### `samtools` (steps 1 and 9)

- `fastq -@ 4`: `-@` sets worker threads for BAM-to-FASTQ conversion.
- `view -b`: `-b` outputs BAM instead of SAM.
- `view -@ <n>` and `sort -@ <n>`: threads for conversion/sorting.
- `sort -o <file>`: output BAM path.
- `index <bam>`: build BAM index (`.bai`).

### `porechop` (step 2)

- `-i <file>`: input FASTQ.
- `-o <file>`: output trimmed FASTQ.
- `-t <n>`: number of threads.

### `NanoFilt` (step 3)

- `-l 100`: keep reads with length >= 100 bp.

### `NanoStat` (steps 4 and 6)

- `--fastq <file>`: input FASTQ file.
- `--name <text>`: sample label used in output.
- `--outdir <dir>`: output directory.
- `--tsv`: write tab-separated output format.

### `kraken2` (steps 5 and 13)

- `--db <dir>`: Kraken2 database directory.
- `--use-names`: print scientific names in classification output.
- `--threads <n>`: parallel threads.
- `--report <file>`: write sample-level clade summary report.
- `--output <file>`: write per-read or per-contig classification output.

### `seqtk` (step 6)

- `subseq <fastq> <id_file>`: extract reads whose IDs are listed in `id_file`.

### `metaMDBG` (step 8)

- `asm`: run assembly mode.
- `--out-dir <dir>`: sample output directory.
- `--in-ont <fastq>`: ONT reads input.
- `--threads <n>`: assembly threads.

### `myloasm` (step 8b)

- `-o <dir>`: output directory.
- `-t <n>`: threads.

### `flye` (step 9)

- `--meta`: metagenome assembly mode.
- `--nano-hq <fastq>`: ONT high-quality reads input.
- `--threads <n>`: threads.
- `-o <dir>`: output directory.

### `minimap2` (step 9)

- `-ax map-ont`: preset for ONT read-to-reference alignment with SAM output.
- `-t <n>`: threads.

### `racon` (step 9)

- `-t <n>`: threads for polishing.

### `seqkit` (steps 10 and 18)

- `fq2fa`: convert FASTQ to FASTA.
- `fq2fa -o <file>`: output FASTA path.
- `stats -a`: include all summary statistics.
- `stats -T`: tabular output.
- `stats -j <n>`: threads.

### `amrfinder` / AMRFinderPlus (step 11)

- `--plus`: include the extended AMRFinderPlus database (AMR + additional marker classes).
- `-n <fasta>`: nucleotide FASTA input.
- `--threads <n>`: threads.
- `-o <file>`: output TSV path.

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

## Final note

If you only change one thing before running: make sure every placeholder path is replaced with real paths for your system. Most failed runs come from path mismatches.

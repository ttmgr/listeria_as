#!/usr/bin/env bash
set -euo pipefail

REMOTE="timthilomaria.reska@hpc-build01"
REMOTE_BASE="/lustre/groups/hpc/urban_lab/projects/tim/foodsafety"
LOCAL_BASE="${1:-./data}"

mkdir -p "$LOCAL_BASE"

FILES=(
  sample_metadata.csv
  processing/listeria/overview/listeria_overview.csv
  processing/read_lengths_filtered_agg_rebuilt.tsv
  processing/read_lengths_filtered_agg.tsv
  processing/amrfinder/overview/amr_reads_overview.csv
  processing/amrfinder/overview/amr_contigs_overview.csv
  processing/stats/assembly_stats_flye.tsv
  processing/stats/assembly_stats_mdbg.tsv
  processing/stats/assembly_stats_myloasm.tsv
  processing/kraken2_csv/reads_classification.csv
  processing/kraken2_csv/flye_contigs_classification.csv
  processing/kraken2_csv/mdbg_contigs_classification.csv
  processing/kraken2_csv/myloasm_contigs_classification.csv
)

printf '%s\n' "${FILES[@]}" | rsync -avz --files-from=- "$REMOTE:$REMOTE_BASE/" "$LOCAL_BASE/"

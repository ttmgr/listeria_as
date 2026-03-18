# Project Overview

## Takeaway

**Tim, for the love of god, do plots locally.**

Running the report builder on your own machine takes seconds, gives you immediate
feedback, and you can iterate on figures without waiting for cluster jobs.

## What this project does

Oxford Nanopore **Adaptive Sampling (AS) vs Normal (N)** sequencing benchmark for
*Listeria monocytogenes* detection on food-contact surfaces. Two experimental
cohorts (**Black** and **Blue**), two sequencing rounds, three extraction methods
(Sponge, Cotton, Zymo), and positive controls (Lm2/Lm4/Lm6).

## What happened on 2026-03-18

### Bugs fixed

1. **Duplicate rows in report tables** — `sample_metadata.csv` has duplicate
   `basename` entries (same barcode appears for 4h and 24h timepoints). The merge
   in `prepare_cohort_tables` created a cartesian product. Fixed by deduplicating
   metadata by `sample` before the merge.

2. **Results mismatch with codex reference** — caused by the same duplicate-row
   issue inflating all metrics. Once deduplication was in place, round-2 reports
   match the earlier codex-only results.

3. **Histogram x-axis ticks** — changed from `100, 1k, 10k, 100k` to
   `1, 10, 100, 1k, 10k, 100k` across all four histogram functions. Bin lower
   bound moved from `log10(50)` to `0` (i.e. 1 bp) so histograms cover the full
   range.

### Features added

4. **Collapsible HTML sections** — every section in the report is now wrapped in
   `<details>/<summary>` so you can collapse sections you're not looking at.

5. **Controls treated separately** — `Control_Lm` samples (Lm2, Lm4, Lm6) are
   excluded from the main cohort analysis (they were inflating/skewing metrics as
   positive controls). They now appear in a dedicated "Controls" section at the
   bottom of each report.

6. **Combined contig taxa CSV** — exported per cohort/round at
   `data/<report>/data/<prefix>_combined_contig_taxa.csv` with columns:
   `cohort, cohort_group, barcode, condition, sample, contig_id, length_bp, taxon, assembler`.

7. **Combined read taxa CSV** — same idea but for per-read kraken2 classifications:
   `data/<report>/data/<prefix>_combined_read_taxa.csv` with columns:
   `cohort, cohort_group, barcode, condition, sample, read_id, length_bp, taxon`.
   Requires `reads_classification.csv` to be downloaded from the cluster first.

8. **Download script updated** — `scripts/download_from_cluster.sh` now also
   fetches `reads_classification.csv`.

## How to regenerate reports

```bash
# 1. Download data from cluster (one-time or when data changes)
bash scripts/download_from_cluster.sh ./data/

# 2. Build all reports (Black + Blue, R1 + R2)
python scripts/build_local_black_report.py ./data/
```

There is only ONE report script. It generates all cohort/round combinations.
Reports land in `data/local_{black,blue}_report_{r1,r2}/report.html`.

## Key files

| File | Purpose |
|------|---------|
| `scripts/build_local_black_report.py` | Main report generator (all cohorts, all rounds) |
| `scripts/download_from_cluster.sh` | rsync wrapper to pull processed data from HPC |
| `data/sample_metadata.csv` | Master sample sheet |
| `data/processing/kraken2_csv/reads_classification.csv` | Per-read kraken2 taxa (~1 GB) |
| `data/processing/kraken2_csv/*_contigs_classification.csv` | Per-contig kraken2 taxa |
| `data/processing/listeria/overview/listeria_overview.csv` | Listeria read/assembly metrics |

#!/usr/bin/env python3
"""
Split downloaded Kraken CSV exports into cohort-specific folders.
Supports round-aware sample names via the 'basename' column in metadata.

Usage:
    python3 split_kraken_csv_by_cohort.py /path/to/downloaded_results/codex

Outputs:
    <base_dir>/processing/kraken2_csv/black/*.csv
    <base_dir>/processing/kraken2_csv/blue/*.csv
"""

from __future__ import annotations

import csv
import sys
from pathlib import Path

import pandas as pd


CSV_NAMES = [
    "reads_classification.csv",
    "flye_contigs_classification.csv",
    "mdbg_contigs_classification.csv",
    "myloasm_contigs_classification.csv",
]

# Chunk size for reading large CSVs (reads_classification.csv can be ~800MB+)
CHUNK_SIZE = 500_000


def cohort_slug(cohort: str) -> str:
    return str(cohort).strip().lower()


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: python3 split_kraken_csv_by_cohort.py <downloaded_base_dir>")
        return 1

    base_dir = Path(sys.argv[1]).resolve()
    meta_path = base_dir / "sample_metadata.csv"
    kraken_dir = base_dir / "processing" / "kraken2_csv"

    meta = pd.read_csv(meta_path)

    # Use 'basename' column if available (round-aware), else fall back to legacy naming
    if "basename" in meta.columns:
        meta["sample"] = meta["basename"]
    else:
        meta["sample"] = meta.apply(lambda r: f"barcode{int(r['barcode']):02d}_{r['condition']}", axis=1)

    cohorts = {
        cohort_slug(cohort): set(sub["sample"].tolist())
        for cohort, sub in meta.groupby("cohort")
        if str(cohort) in {"Black", "Blue"}
    }

    # Also build round-specific cohorts for per-round splitting
    rounds = sorted(meta["round"].dropna().unique()) if "round" in meta.columns else []
    round_cohorts = {}
    if len(rounds) > 1:
        for rnd in rounds:
            rnd_int = int(rnd)
            rnd_meta = meta[meta["round"] == rnd]
            for cohort, sub in rnd_meta.groupby("cohort"):
                if str(cohort) not in {"Black", "Blue"}:
                    continue
                slug = cohort_slug(cohort)
                key = f"r{rnd_int}_{slug}"
                round_cohorts[key] = set(sub["sample"].tolist())

    all_groups = {**cohorts, **round_cohorts}
    for slug in all_groups:
        (kraken_dir / slug).mkdir(parents=True, exist_ok=True)

    for csv_name in CSV_NAMES:
        src = kraken_dir / csv_name
        if not src.exists():
            continue

        # Use chunked reading for large files
        is_large = csv_name == "reads_classification.csv"

        if is_large:
            # Initialize output files with headers
            first_chunk = True
            writers = {}
            for chunk in pd.read_csv(src, dtype={"sample": str}, quoting=csv.QUOTE_MINIMAL, chunksize=CHUNK_SIZE):
                if "sample" not in chunk.columns:
                    print(f"Skipping {csv_name}: no sample column")
                    break
                for slug, samples in all_groups.items():
                    out = kraken_dir / slug / csv_name
                    filtered = chunk[chunk["sample"].isin(samples)]
                    if len(filtered) > 0:
                        filtered.to_csv(out, index=False, mode="a" if not first_chunk else "w",
                                        header=first_chunk)
                first_chunk = False
            for slug in all_groups:
                out = kraken_dir / slug / csv_name
                if out.exists():
                    print(f"Wrote {out}")
        else:
            df = pd.read_csv(src, dtype={"sample": str}, quoting=csv.QUOTE_MINIMAL)
            if "sample" not in df.columns:
                print(f"Skipping {csv_name}: no sample column")
                continue
            for slug, samples in all_groups.items():
                out = kraken_dir / slug / csv_name
                df[df["sample"].isin(samples)].to_csv(out, index=False)
                print(f"Wrote {out}")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
Split downloaded AMR overview CSVs into cohort-specific folders.

Usage:
    python3 split_amr_overview_by_cohort.py /path/to/downloaded_results/codex

Outputs:
    <base_dir>/processing/amrfinder/overview/black/amr_reads_overview.csv
    <base_dir>/processing/amrfinder/overview/black/amr_contigs_overview.csv
    <base_dir>/processing/amrfinder/overview/blue/amr_reads_overview.csv
    <base_dir>/processing/amrfinder/overview/blue/amr_contigs_overview.csv
"""

from __future__ import annotations

import re
import sys
from pathlib import Path

import pandas as pd


CSV_NAMES = ["amr_reads_overview.csv", "amr_contigs_overview.csv"]


def cohort_slug(cohort: str) -> str:
    return str(cohort).strip().lower()


def barcode_to_int(value: object) -> int | None:
    match = re.search(r"(\d+)", str(value))
    return int(match.group(1)) if match else None


def extract_round(sample_name: str) -> int | None:
    """Extract round number from sample name like r1_barcode03_AS -> 1."""
    m = re.match(r"r(\d+)_", str(sample_name))
    return int(m.group(1)) if m else None


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: python3 split_amr_overview_by_cohort.py <downloaded_base_dir>")
        return 1

    base_dir = Path(sys.argv[1]).resolve()
    meta = pd.read_csv(base_dir / "sample_metadata.csv")
    overview_dir = base_dir / "processing" / "amrfinder" / "overview"

    # Use basename column if available (round-aware)
    if "basename" in meta.columns:
        meta["_sample"] = meta["basename"]
    else:
        meta["_sample"] = meta.apply(lambda r: f"barcode{int(r['barcode']):02d}_{r['condition']}", axis=1)

    # Build cohort -> set of basenames
    cohort_samples = {}
    for cohort, sub in meta.groupby("cohort"):
        if str(cohort) not in {"Black", "Blue"}:
            continue
        cohort_samples[cohort_slug(cohort)] = set(sub["_sample"].tolist())

    # Also build barcode-based lookup for legacy files without Sample column
    cohort_barcodes = {}
    for cohort, sub in meta.groupby("cohort"):
        if str(cohort) not in {"Black", "Blue"}:
            continue
        cohort_barcodes[cohort_slug(cohort)] = set(sub["barcode"].dropna().astype(int).tolist())

    # Build round-specific groups
    rounds = sorted(meta["round"].dropna().unique()) if "round" in meta.columns else []
    round_cohort_samples = {}
    round_cohort_barcodes = {}
    if len(rounds) > 1:
        for rnd in rounds:
            rnd_meta = meta[meta["round"] == rnd]
            for cohort, sub in rnd_meta.groupby("cohort"):
                if str(cohort) not in {"Black", "Blue"}:
                    continue
                slug = cohort_slug(cohort)
                key = f"r{int(rnd)}_{slug}"
                round_cohort_samples[key] = set(sub["_sample"].tolist())
                round_cohort_barcodes[key] = set(sub["barcode"].dropna().astype(int).tolist())

    all_sample_groups = {**cohort_samples, **round_cohort_samples}
    all_barcode_groups = {**cohort_barcodes, **round_cohort_barcodes}

    for slug in all_sample_groups:
        (overview_dir / slug).mkdir(parents=True, exist_ok=True)

    for csv_name in CSV_NAMES:
        src = overview_dir / csv_name
        if not src.exists():
            continue
        df = pd.read_csv(src)

        # Try to match by Sample column first, fall back to Barcode
        if "Sample" in df.columns:
            for slug, samples in all_sample_groups.items():
                out = overview_dir / slug / csv_name
                df[df["Sample"].isin(samples)].to_csv(out, index=False)
                print(f"Wrote {out}")
        elif "Barcode" in df.columns:
            df["_barcode_num"] = df["Barcode"].map(barcode_to_int)
            for slug, barcodes in all_barcode_groups.items():
                out = overview_dir / slug / csv_name
                df[df["_barcode_num"].isin(barcodes)].drop(columns=["_barcode_num"]).to_csv(out, index=False)
                print(f"Wrote {out}")
        else:
            print(f"Skipping {csv_name}: no Sample or Barcode column")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

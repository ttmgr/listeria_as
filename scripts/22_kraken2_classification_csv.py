#!/usr/bin/env python3
"""
Step 22: Kraken2 Classification CSV Export.
Purpose: parse Kraken2 per-read and per-contig classification files and produce
         clean CSV tables with columns: id, length_bp, taxon.

Produces 4 CSV files (one for reads, one per assembler):
  - reads_classification.csv
  - flye_contigs_classification.csv
  - mdbg_contigs_classification.csv
  - myloasm_contigs_classification.csv

Each row corresponds to one read or contig with the Kraken2 taxon assignment.

Usage:
    python3 22_kraken2_classification_csv.py <base_dir>

    base_dir: project root that contains processing/kraken2/ and processing/kraken2_contigs/
"""

from __future__ import annotations

import csv
import os
import re
import sys
from pathlib import Path
from typing import List, Optional, Tuple


# ---------------------------------------------------------------------------
# Kraken2 classified output format (tab-separated):
#   col 0: C/U  (classified / unclassified)
#   col 1: read/contig ID
#   col 2: taxon name (with --use-names) or taxid
#   col 3: sequence length  (e.g. "1234" or "100|100" for paired)
#   col 4: kmer mapping info
# ---------------------------------------------------------------------------


def parse_length(raw: str) -> int:
    """Parse Kraken2 length field — handles both single int and paired '100|100'."""
    nums = re.findall(r"\d+", raw)
    if not nums:
        return 0
    return sum(int(n) for n in nums)


def parse_taxon(raw: str) -> str:
    """Clean up the taxon field — strip trailing '(taxid NNN)' if present."""
    raw = raw.strip()
    m = re.match(r"^(.*?)\s*\(taxid\s+\d+\)\s*$", raw)
    if m:
        return m.group(1).strip()
    return raw


def derive_sample_from_filename(filename: str) -> str:
    """Extract sample name from classified file name.

    Examples:
        classified_r1_barcode03_AS.txt  -> r1_barcode03_AS
        classified_barcode03_AS.txt     -> barcode03_AS
        classified_flye_r1_barcode03_AS.txt -> r1_barcode03_AS
        classified_mdbg_barcode03_AS.txt    -> barcode03_AS
    """
    stem = filename.replace(".txt", "")
    # Remove prefixes: classified_, classified_flye_, classified_mdbg_, classified_myloasm_
    for prefix in ["classified_myloasm_", "classified_mdbg_", "classified_flye_", "classified_"]:
        if stem.startswith(prefix):
            stem = stem[len(prefix):]
            break
    return stem


def parse_classified_file(path: Path) -> List[Tuple[str, int, str, str]]:
    """Parse one Kraken2 classified_*.txt file.

    Returns list of (id, length_bp, taxon, sample) tuples.
    """
    sample = derive_sample_from_filename(path.name)
    rows: List[Tuple[str, int, str, str]] = []
    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            parts = line.split("\t", 4)
            if len(parts) < 4:
                continue
            status = parts[0].strip()
            seq_id = parts[1].strip()
            taxon_raw = parts[2]
            length_raw = parts[3]

            taxon = parse_taxon(taxon_raw) if status == "C" else "unclassified"
            length = parse_length(length_raw)
            rows.append((seq_id, length, taxon, sample))
    return rows


def write_csv(rows: List[Tuple[str, int, str, str]], out_path: Path, id_label: str) -> int:
    """Write rows to a CSV with header id_label, length_bp, taxon, sample."""
    out_path.parent.mkdir(parents=True, exist_ok=True)
    with out_path.open("w", newline="", encoding="utf-8") as fh:
        writer = csv.writer(fh)
        writer.writerow([id_label, "length_bp", "taxon", "sample"])
        writer.writerows(rows)
    return len(rows)


def collect_files(directory: Path, pattern: str) -> List[Path]:
    """Recursively find files matching a glob pattern."""
    if not directory.exists():
        return []
    return sorted(directory.rglob(pattern))


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------

def main() -> int:
    if len(sys.argv) < 2:
        print("Usage: python3 22_kraken2_classification_csv.py <base_dir>")
        return 1

    base_dir = Path(sys.argv[1]).resolve()
    out_dir = base_dir / "processing" / "kraken2_csv"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Define sources
    sources = {
        "reads": {
            "search_dir": base_dir / "processing" / "kraken2",
            "pattern": "classified_*.txt",
            "out_file": "reads_classification.csv",
            "id_label": "read_id",
        },
        "flye": {
            "search_dir": base_dir / "processing" / "kraken2_contigs" / "flye",
            "pattern": "classified_flye_*.txt",
            "out_file": "flye_contigs_classification.csv",
            "id_label": "contig_id",
        },
        "mdbg": {
            "search_dir": base_dir / "processing" / "kraken2_contigs" / "mdbg",
            "pattern": "classified_mdbg_*.txt",
            "out_file": "mdbg_contigs_classification.csv",
            "id_label": "contig_id",
        },
        "myloasm": {
            "search_dir": base_dir / "processing" / "kraken2_contigs" / "myloasm",
            "pattern": "classified_myloasm_*.txt",
            "out_file": "myloasm_contigs_classification.csv",
            "id_label": "contig_id",
        },
    }

    total_exported = 0

    for name, cfg in sources.items():
        search_dir = cfg["search_dir"]
        files = collect_files(search_dir, cfg["pattern"])
        print(f"[{name}] Found {len(files)} classified files in {search_dir}")

        all_rows: List[Tuple[str, int, str, str]] = []
        for f in files:
            rows = parse_classified_file(f)
            all_rows.extend(rows)
            print(f"  Parsed {f.name}: {len(rows)} records")

        if all_rows:
            out_path = out_dir / cfg["out_file"]
            n = write_csv(all_rows, out_path, cfg["id_label"])
            print(f"  Saved: {out_path} ({n} rows)")
            total_exported += n
        else:
            print(f"  No records found for {name}, skipping CSV.")

    # Summary
    print()
    print("=" * 60)
    print("  KRAKEN2 CLASSIFICATION CSV EXPORT SUMMARY")
    print("=" * 60)
    print(f"  Output directory: {out_dir}")
    print(f"  Total records exported: {total_exported:,}")
    for name, cfg in sources.items():
        p = out_dir / cfg["out_file"]
        if p.exists():
            # Count lines (minus header)
            with p.open() as fh:
                n = sum(1 for _ in fh) - 1
            print(f"  {cfg['out_file']}: {n:,} rows")
        else:
            print(f"  {cfg['out_file']}: not created")
    print("=" * 60)

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

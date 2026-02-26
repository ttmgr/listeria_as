#!/usr/bin/env python3
"""
Export Kraken2 read/contig classifications to spreadsheet-ready tables.

Input files:
  - classified_*.txt
  - report_*.txt

Works for:
  - Reads Kraken2 outputs (typically under kraken2/)
  - Contig Kraken2 outputs (typically under kraken2_contigs/ and flye/mdbg subdirs)

Outputs:
  - CSV tables in <output_dir>
  - XLSX files for each CSV (using pandas Excel writer or ssconvert fallback)
  - Combined workbook all_kraken2_tables.xlsx (best effort)
  - export_manifest.csv

Usage:
  python3 21_kraken2_to_spreadsheets.py /path/to/base_dir
"""

from __future__ import annotations

import argparse
import re
import shutil
import subprocess
from collections import defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple

import pandas as pd


CLASSIFIED_RE = re.compile(r"^(classified)_(?:(flye|mdbg|metaMDBG)_)?(.+)\.txt$", re.IGNORECASE)
REPORT_RE = re.compile(r"^(report)_(?:(flye|mdbg|metaMDBG)_)?(.+)\.txt$", re.IGNORECASE)
TAXID_RE = re.compile(r"^(.*)\s+\(taxid\s+(\d+)\)\s*$")


def normalize_assembler(val: Optional[str]) -> str:
    if not val:
        return ""
    low = val.lower()
    if low == "metamdbg":
        return "mdbg"
    return low


def find_files(root: Path, pattern: str) -> List[Path]:
    if not root.exists():
        return []
    return sorted(p for p in root.rglob(pattern) if p.is_file())


def infer_metadata(path: Path, source: str) -> Tuple[str, str]:
    name = path.name
    assembler = ""
    sample = path.stem

    m = CLASSIFIED_RE.match(name) or REPORT_RE.match(name)
    if m:
        assembler = normalize_assembler(m.group(2))
        sample = m.group(3)

    # Fallback from directory names for contigs.
    if source == "contigs" and not assembler:
        parent = str(path.parent).lower()
        if "/flye/" in parent or parent.endswith("/flye"):
            assembler = "flye"
        elif "mdbg" in parent:
            assembler = "mdbg"

    if source == "reads":
        assembler = "reads"
    elif not assembler:
        assembler = "unknown"

    return sample, assembler


def parse_taxon_field(taxon_raw: str) -> Tuple[str, Optional[int]]:
    taxon_raw = taxon_raw.strip()
    m = TAXID_RE.match(taxon_raw)
    if m:
        return m.group(1).strip(), int(m.group(2))
    if taxon_raw.isdigit():
        return taxon_raw, int(taxon_raw)
    return taxon_raw, None


def parse_length_field(length_raw: str) -> Optional[int]:
    nums = re.findall(r"\d+", length_raw)
    if not nums:
        return None
    # Kraken2 single-end outputs one value. For paired formats like "100|100",
    # summing keeps a consistent total-bases interpretation.
    return sum(int(x) for x in nums)


def parse_report_line(line: str) -> Optional[Dict[str, object]]:
    stripped = line.rstrip("\n")
    if not stripped.strip():
        return None

    parts = stripped.split("\t")
    if len(parts) >= 6:
        pct, clade_cnt, taxon_cnt, rank_code, taxid, sci = parts[:6]
    else:
        parts = stripped.split(None, 5)
        if len(parts) < 6:
            return None
        pct, clade_cnt, taxon_cnt, rank_code, taxid, sci = parts

    sci_raw = sci.rstrip("\n")
    indent = len(sci_raw) - len(sci_raw.lstrip(" "))
    sci_clean = sci_raw.strip()

    def to_int(val: str) -> Optional[int]:
        try:
            return int(float(val))
        except Exception:
            return None

    def to_float(val: str) -> Optional[float]:
        try:
            return float(val)
        except Exception:
            return None

    return {
        "percentage": to_float(pct),
        "clade_count": to_int(clade_cnt),
        "taxon_count": to_int(taxon_cnt),
        "rank_code": rank_code.strip(),
        "taxid": to_int(taxid),
        "scientific_name": sci_clean,
        "indent_spaces": indent,
    }


def process_classified_file(path: Path, source: str) -> Tuple[Dict[str, object], List[Dict[str, object]]]:
    sample, assembler = infer_metadata(path, source)
    total = 0
    classified = 0
    unclassified = 0
    bad_lines = 0

    len_sum_all = 0
    len_n_all = 0
    len_sum_c = 0
    len_n_c = 0
    len_sum_u = 0
    len_n_u = 0

    tax_counts: Dict[Tuple[str, Optional[int], str], int] = defaultdict(int)

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            if not line.strip():
                continue
            parts = line.rstrip("\n").split("\t", 4)
            if len(parts) < 4:
                bad_lines += 1
                continue

            status = parts[0].strip()
            taxon_raw = parts[2]
            length_raw = parts[3]
            tax_name, taxid = parse_taxon_field(taxon_raw)
            seq_len = parse_length_field(length_raw)

            total += 1
            tax_counts[(status, taxid, tax_name)] += 1

            if status == "C":
                classified += 1
            elif status == "U":
                unclassified += 1

            if seq_len is not None:
                len_sum_all += seq_len
                len_n_all += 1
                if status == "C":
                    len_sum_c += seq_len
                    len_n_c += 1
                elif status == "U":
                    len_sum_u += seq_len
                    len_n_u += 1

    def mean(sum_v: int, n_v: int) -> Optional[float]:
        if n_v == 0:
            return None
        return round(sum_v / n_v, 3)

    summary = {
        "source": source,
        "assembler": assembler,
        "sample": sample,
        "file_name": path.name,
        "file_path": str(path),
        "total_records": total,
        "classified_records": classified,
        "unclassified_records": unclassified,
        "classified_pct": round((classified / total) * 100, 4) if total else 0.0,
        "unclassified_pct": round((unclassified / total) * 100, 4) if total else 0.0,
        "mean_length_bp": mean(len_sum_all, len_n_all),
        "mean_length_classified_bp": mean(len_sum_c, len_n_c),
        "mean_length_unclassified_bp": mean(len_sum_u, len_n_u),
        "bad_lines": bad_lines,
    }

    tax_rows = []
    for (status, taxid, tax_name), cnt in tax_counts.items():
        tax_rows.append(
            {
                "source": source,
                "assembler": assembler,
                "sample": sample,
                "status": status,
                "taxid": taxid,
                "taxon_name": tax_name,
                "count": cnt,
                "pct_of_sample": round((cnt / total) * 100, 4) if total else 0.0,
                "file_name": path.name,
            }
        )

    return summary, tax_rows


def process_report_file(path: Path, source: str) -> List[Dict[str, object]]:
    sample, assembler = infer_metadata(path, source)
    rows: List[Dict[str, object]] = []
    bad_lines = 0

    with path.open("r", encoding="utf-8", errors="replace") as fh:
        for line in fh:
            parsed = parse_report_line(line)
            if parsed is None:
                if line.strip():
                    bad_lines += 1
                continue
            rows.append(
                {
                    "source": source,
                    "assembler": assembler,
                    "sample": sample,
                    "file_name": path.name,
                    "file_path": str(path),
                    **parsed,
                }
            )

    if bad_lines:
        rows.append(
            {
                "source": source,
                "assembler": assembler,
                "sample": sample,
                "file_name": path.name,
                "file_path": str(path),
                "percentage": None,
                "clade_count": None,
                "taxon_count": None,
                "rank_code": "_META",
                "taxid": None,
                "scientific_name": f"BAD_LINES={bad_lines}",
                "indent_spaces": None,
            }
        )
    return rows


def write_csv(df: pd.DataFrame, out_path: Path) -> Path:
    out_path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(out_path, index=False)
    return out_path


def csv_to_xlsx(csv_path: Path, xlsx_path: Path) -> bool:
    # Try pandas Excel engines first.
    for engine in ("openpyxl", "xlsxwriter"):
        try:
            df = pd.read_csv(csv_path)
            df.to_excel(xlsx_path, index=False, engine=engine)
            return True
        except Exception:
            continue

    # Fallback to ssconvert if available.
    if shutil.which("ssconvert"):
        try:
            subprocess.run(
                ["ssconvert", str(csv_path), str(xlsx_path)],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            return True
        except Exception:
            return False
    return False


def merge_xlsx_files(xlsx_files: Iterable[Path], out_path: Path) -> bool:
    xlsx_files = list(xlsx_files)
    if not xlsx_files:
        return False

    if shutil.which("ssconvert"):
        try:
            subprocess.run(
                ["ssconvert", "--merge-to", str(out_path)] + [str(p) for p in xlsx_files],
                check=True,
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
            )
            return True
        except Exception:
            return False
    return False


def make_top_taxa(df_taxon: pd.DataFrame, top_n: int = 20) -> pd.DataFrame:
    if df_taxon.empty:
        return df_taxon
    cdf = df_taxon[df_taxon["status"] == "C"].copy()
    if cdf.empty:
        return cdf
    cdf = cdf.sort_values(["source", "assembler", "sample", "count"], ascending=[True, True, True, False])
    cdf["rank_in_sample"] = cdf.groupby(["source", "assembler", "sample"])["count"].rank(
        method="first", ascending=False
    )
    cdf = cdf[cdf["rank_in_sample"] <= top_n].copy()
    cdf["rank_in_sample"] = cdf["rank_in_sample"].astype(int)
    return cdf


def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(description="Convert Kraken2 outputs into spreadsheet-ready tables.")
    p.add_argument("base_dir", nargs="?", default=".", help="Base directory that contains kraken2/ and kraken2_contigs/")
    p.add_argument("--reads-dir", default="kraken2", help="Reads Kraken2 directory (default: kraken2)")
    p.add_argument("--contigs-dir", default="kraken2_contigs", help="Contigs Kraken2 directory (default: kraken2_contigs)")
    p.add_argument(
        "--output-dir",
        default="processing/kraken2_spreadsheets",
        help="Output directory (default: processing/kraken2_spreadsheets)",
    )
    return p


def main() -> int:
    args = build_parser().parse_args()
    base_dir = Path(args.base_dir).resolve()
    reads_dir = (base_dir / args.reads_dir).resolve()
    contigs_dir = (base_dir / args.contigs_dir).resolve()
    out_dir = (base_dir / args.output_dir).resolve()
    out_dir.mkdir(parents=True, exist_ok=True)

    reads_classified = find_files(reads_dir, "classified_*.txt")
    reads_report = find_files(reads_dir, "report_*.txt")
    contigs_classified = find_files(contigs_dir, "classified_*.txt")
    contigs_report = find_files(contigs_dir, "report_*.txt")

    print(f"Base dir: {base_dir}")
    print(f"Reads dir: {reads_dir}")
    print(f"Contigs dir: {contigs_dir}")
    print(f"Found reads classified: {len(reads_classified)}")
    print(f"Found reads report: {len(reads_report)}")
    print(f"Found contigs classified: {len(contigs_classified)}")
    print(f"Found contigs report: {len(contigs_report)}")

    if not (reads_classified or reads_report or contigs_classified or contigs_report):
        print("No Kraken2 files found. Nothing to export.")
        return 0

    classified_summaries: List[Dict[str, object]] = []
    taxon_rows: List[Dict[str, object]] = []
    report_rows: List[Dict[str, object]] = []

    for f in reads_classified:
        s, t = process_classified_file(f, "reads")
        classified_summaries.append(s)
        taxon_rows.extend(t)
        print(f"Parsed classified reads: {f.name}")

    for f in contigs_classified:
        s, t = process_classified_file(f, "contigs")
        classified_summaries.append(s)
        taxon_rows.extend(t)
        print(f"Parsed classified contigs: {f.name}")

    for f in reads_report:
        rows = process_report_file(f, "reads")
        report_rows.extend(rows)
        print(f"Parsed report reads: {f.name}")

    for f in contigs_report:
        rows = process_report_file(f, "contigs")
        report_rows.extend(rows)
        print(f"Parsed report contigs: {f.name}")

    csv_paths: List[Path] = []

    if classified_summaries:
        df_cls = pd.DataFrame(classified_summaries).sort_values(["source", "assembler", "sample"])
        csv_paths.append(write_csv(df_cls, out_dir / "classified_sample_summary.csv"))
    else:
        df_cls = pd.DataFrame()

    if taxon_rows:
        df_tax = pd.DataFrame(taxon_rows).sort_values(
            ["source", "assembler", "sample", "status", "count"],
            ascending=[True, True, True, True, False],
        )
        csv_paths.append(write_csv(df_tax, out_dir / "classified_taxon_counts.csv"))

        df_top = make_top_taxa(df_tax, top_n=20)
        csv_paths.append(write_csv(df_top, out_dir / "classified_top20_taxa_per_sample.csv"))
    else:
        df_tax = pd.DataFrame()

    if report_rows:
        df_rep = pd.DataFrame(report_rows).sort_values(
            ["source", "assembler", "sample", "rank_code", "taxon_count"],
            ascending=[True, True, True, True, False],
        )
        csv_paths.append(write_csv(df_rep, out_dir / "report_entries.csv"))

        species = df_rep[df_rep["rank_code"] == "S"].copy()
        if not species.empty:
            species = species.sort_values(
                ["source", "assembler", "sample", "taxon_count"],
                ascending=[True, True, True, False],
            )
            species["rank_in_sample"] = species.groupby(["source", "assembler", "sample"])["taxon_count"].rank(
                method="first", ascending=False
            )
            species = species[species["rank_in_sample"] <= 20].copy()
            species["rank_in_sample"] = species["rank_in_sample"].astype(int)
            csv_paths.append(write_csv(species, out_dir / "report_top20_species_per_sample.csv"))
    else:
        df_rep = pd.DataFrame()

    manifest_rows = []
    xlsx_paths: List[Path] = []

    for csv_path in csv_paths:
        xlsx_path = csv_path.with_suffix(".xlsx")
        ok = csv_to_xlsx(csv_path, xlsx_path)
        manifest_rows.append(
            {
                "csv_file": csv_path.name,
                "xlsx_file": xlsx_path.name if ok else "",
                "xlsx_created": ok,
                "rows": len(pd.read_csv(csv_path)),
            }
        )
        if ok:
            xlsx_paths.append(xlsx_path)
        print(f"Exported CSV: {csv_path.name}")
        print(f"Exported XLSX: {xlsx_path.name if ok else 'FAILED'}")

    combined_xlsx = out_dir / "all_kraken2_tables.xlsx"
    combined_ok = merge_xlsx_files(xlsx_paths, combined_xlsx)
    print(f"Combined workbook: {'created' if combined_ok else 'not created'}")

    manifest = pd.DataFrame(manifest_rows)
    if combined_ok:
        manifest = pd.concat(
            [
                manifest,
                pd.DataFrame(
                    [
                        {
                            "csv_file": "",
                            "xlsx_file": combined_xlsx.name,
                            "xlsx_created": True,
                            "rows": "",
                        }
                    ]
                ),
            ],
            ignore_index=True,
        )
    write_csv(manifest, out_dir / "export_manifest.csv")

    print(f"All Kraken2 tables exported to: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

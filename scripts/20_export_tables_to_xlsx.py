#!/usr/bin/env python3
"""
Step 20: Export pipeline tables as spreadsheets.

Exports every CSV/TSV table under <base_dir>/processing to:
  - Individual XLSX files in <base_dir>/processing/spreadsheets/
  - One combined workbook: <base_dir>/processing/spreadsheets/all_tables.xlsx
  - Manifest: <base_dir>/processing/spreadsheets/spreadsheet_manifest.csv

Usage:
  python3 20_export_tables_to_xlsx.py <base_dir>
"""

import csv
import re
import shutil
import subprocess
import sys
from pathlib import Path

import pandas as pd


def to_safe_filename(rel_path: Path) -> str:
    stem = str(rel_path.with_suffix("")).replace("/", "__")
    return re.sub(r"[^A-Za-z0-9._-]+", "_", stem) + ".xlsx"


def to_sheet_name(rel_path: Path, used: set[str]) -> str:
    base = re.sub(r"[^A-Za-z0-9_]+", "_", str(rel_path.with_suffix("")))
    base = base[:31] if base else "sheet"
    name = base
    i = 1
    while name in used:
        suffix = f"_{i}"
        name = f"{base[:31 - len(suffix)]}{suffix}"
        i += 1
    used.add(name)
    return name


def normalize_header(raw_header: list[str], width: int) -> list[str]:
    cols = list(raw_header[:width])
    if len(cols) < width:
        cols.extend([""] * (width - len(cols)))

    used: dict[str, int] = {}
    out = []
    for i, col in enumerate(cols):
        name = (col or "").strip() or f"col_{i + 1}"
        if name in used:
            used[name] += 1
            name = f"{name}_{used[name]}"
        else:
            used[name] = 0
        out.append(name)
    return out


def load_table(path: Path) -> pd.DataFrame:
    delimiter = "\t" if path.suffix.lower() == ".tsv" else ","
    with path.open(newline="", encoding="utf-8", errors="replace") as f:
        rows = list(csv.reader(f, delimiter=delimiter))

    if not rows:
        return pd.DataFrame()

    width = max(len(r) for r in rows)
    rows = [r + [""] * (width - len(r)) for r in rows]
    header = normalize_header(rows[0], width)
    return pd.DataFrame(rows[1:], columns=header)


def write_xlsx_pandas(df: pd.DataFrame, out_xlsx: Path) -> bool:
    for engine in ("openpyxl", "xlsxwriter"):
        try:
            df.to_excel(out_xlsx, index=False, engine=engine)
            return True
        except Exception:
            continue
    return False


def write_xlsx_ssconvert(src_table: Path, out_xlsx: Path) -> bool:
    if shutil.which("ssconvert") is None:
        return False
    try:
        subprocess.run(
            ["ssconvert", str(src_table), str(out_xlsx)],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return True
    except Exception:
        return False


def create_combined_workbook_pandas(tables: list[Path], processing_dir: Path, combined_path: Path) -> bool:
    try:
        used_sheet_names: set[str] = set()
        with pd.ExcelWriter(combined_path, engine="openpyxl") as writer:
            for table in tables:
                rel = table.relative_to(processing_dir)
                df = load_table(table)
                sheet = to_sheet_name(rel, used_sheet_names)
                df.to_excel(writer, sheet_name=sheet, index=False)
        return True
    except Exception:
        return False


def create_combined_workbook_ssconvert(individual_xlsx: list[Path], combined_path: Path) -> bool:
    if shutil.which("ssconvert") is None:
        return False
    if not individual_xlsx:
        return False
    try:
        subprocess.run(
            ["ssconvert", "--merge-to", str(combined_path)] + [str(p) for p in individual_xlsx],
            check=True,
            stdout=subprocess.DEVNULL,
            stderr=subprocess.DEVNULL,
        )
        return True
    except Exception:
        return False


def main() -> int:
    if len(sys.argv) < 2:
        print("Usage: python3 20_export_tables_to_xlsx.py <base_dir>")
        return 1

    base_dir = Path(sys.argv[1]).resolve()
    processing_dir = base_dir / "processing"
    out_dir = processing_dir / "spreadsheets"
    out_dir.mkdir(parents=True, exist_ok=True)

    tables = sorted(
        p for p in processing_dir.rglob("*")
        if p.is_file() and p.suffix.lower() in {".csv", ".tsv"} and out_dir not in p.parents
    )

    if not tables:
        print(f"No CSV/TSV tables found under: {processing_dir}")
        return 0

    manifest_rows = []
    individual_xlsx = []
    used_ssconvert_fallback = False

    # Export individual spreadsheets.
    for table in tables:
        rel = table.relative_to(processing_dir)
        df = load_table(table)
        out_xlsx = out_dir / to_safe_filename(rel)

        wrote = write_xlsx_pandas(df, out_xlsx)
        if not wrote:
            wrote = write_xlsx_ssconvert(table, out_xlsx)
            if wrote:
                used_ssconvert_fallback = True

        if not wrote:
            print(f"FAILED: {rel} -> could not write XLSX (no supported writer found)")
            continue

        individual_xlsx.append(out_xlsx)
        manifest_rows.append(
            {
                "source_table": str(rel),
                "rows": len(df),
                "columns": len(df.columns),
                "spreadsheet_file": out_xlsx.name,
            }
        )
        print(f"Exported: {rel} -> {out_xlsx.name}")

    # Export one combined workbook.
    combined_path = out_dir / "all_tables.xlsx"
    combined_ok = create_combined_workbook_pandas(tables, processing_dir, combined_path)
    if not combined_ok:
        combined_ok = create_combined_workbook_ssconvert(individual_xlsx, combined_path)
        if combined_ok:
            used_ssconvert_fallback = True

    if combined_ok:
        print(f"Exported combined workbook: {combined_path.name}")
    else:
        print("WARNING: Combined workbook could not be created.")

    manifest = pd.DataFrame(manifest_rows).sort_values("source_table")
    manifest_path = out_dir / "spreadsheet_manifest.csv"
    manifest.to_csv(manifest_path, index=False)
    print(f"Saved manifest: {manifest_path.name}")

    if used_ssconvert_fallback:
        print("Note: Used ssconvert fallback for XLSX export.")

    print(f"All spreadsheets saved to: {out_dir}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

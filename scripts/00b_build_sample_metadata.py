#!/usr/bin/env python3
"""
Step 0b: Build a consolidated sample_metadata.csv from the master sample sheet.

Reads samplesheets/sample_sheet_master.csv and produces a flat metadata table
with one row per BAM file (i.e., two rows per sample -- one AS, one N).

Column names are compatible with both the cluster report (17_generate_report_v2.py)
and the local cohort reports (build_local_black_report.py).

Usage:
    python3 scripts/00b_build_sample_metadata.py /path/to/foodsafety

Output: /path/to/foodsafety/sample_metadata.csv
"""
import csv
import os
import sys


def main():
    if len(sys.argv) < 2:
        print("Usage: python3 00b_build_sample_metadata.py <work_dir>")
        sys.exit(1)

    work_dir = sys.argv[1]
    master_csv = os.path.join(work_dir, "samplesheets", "sample_sheet_master.csv")
    output_csv = os.path.join(work_dir, "sample_metadata.csv")

    if not os.path.isfile(master_csv):
        print(f"ERROR: Master sample sheet not found: {master_csv}")
        sys.exit(1)

    out_header = [
        "sample_id", "original_sample_id", "round", "barcode", "barcode_label",
        "condition", "cohort", "group", "swab_type", "kit",
        "dna_concentration_ng_ul", "bam_path", "basename", "comment",
    ]

    rows = []
    with open(master_csv, newline="", encoding="utf-8-sig") as f:
        reader = csv.DictReader(f)
        for row in reader:
            round_num = row.get("Round", "").strip()
            sample_id = row.get("Sample_ID", "").strip()
            barcode = row.get("Barcode", "").strip()
            label_colour = row.get("Label_colour", "").strip()
            enrichment = row.get("Enrichment_condition", "").strip()
            swab_type = row.get("Swab_type", "").strip()
            kit_used = row.get("Kit_used", "").strip()
            quanti = row.get("QuantiFluor_ng_ul", "").strip()
            comment = row.get("Comment", "").strip()
            bam_as = row.get("BAM_AS", "").strip()
            bam_n = row.get("BAM_N", "").strip()

            # Derive cohort from label_colour (Black, Blue, Control, etc.)
            cohort = label_colour if label_colour else "Unknown"

            # Derive group: e.g. "Black_A" based on swab/kit mapping
            # Method letter: Sponge→A, Cotton→C, Zymo→D, else use swab initial
            method_letter = {
                "Sponge": "A", "Cotton": "C", "Zymo": "D",
            }.get(swab_type, swab_type[:1].upper() if swab_type else "X")

            # Control samples get their own group
            if cohort == "Control" or "Lm" in sample_id:
                group = f"Control_Lm"
            else:
                group = f"{cohort}_{method_letter}"

            barcode_label = f"barcode{int(barcode):02d}" if barcode.isdigit() else barcode

            def map_bam_path(bam_col_value):
                if not bam_col_value:
                    return ""
                parts = bam_col_value.split("/", 1)
                if len(parts) == 2:
                    dir_part, file_part = parts
                    rnum = dir_part.replace("round", "").replace("_bam", "")
                    return f"listeria_{rnum}/{file_part}"
                return bam_col_value

            def derive_basename(bam_path):
                """Same logic as the shell derive_basename."""
                if not bam_path:
                    return ""
                dir_name = os.path.basename(os.path.dirname(bam_path))
                file_name = os.path.splitext(os.path.basename(bam_path))[0]
                round_n = dir_name.replace("listeria_", "")
                return f"r{round_n}_{file_name}"

            # Emit one row for AS
            if bam_as:
                mapped_as = map_bam_path(bam_as)
                rows.append({
                    "sample_id": sample_id,
                    "original_sample_id": sample_id,
                    "round": round_num,
                    "barcode": barcode,
                    "barcode_label": barcode_label,
                    "condition": "AS",
                    "cohort": cohort,
                    "group": group,
                    "swab_type": swab_type,
                    "kit": kit_used,
                    "dna_concentration_ng_ul": quanti,
                    "bam_path": mapped_as,
                    "basename": derive_basename(mapped_as),
                    "comment": comment,
                })

            # Emit one row for N
            if bam_n:
                mapped_n = map_bam_path(bam_n)
                rows.append({
                    "sample_id": sample_id,
                    "original_sample_id": sample_id,
                    "round": round_num,
                    "barcode": barcode,
                    "barcode_label": barcode_label,
                    "condition": "N",
                    "cohort": cohort,
                    "group": group,
                    "swab_type": swab_type,
                    "kit": kit_used,
                    "dna_concentration_ng_ul": quanti,
                    "bam_path": mapped_n,
                    "basename": derive_basename(mapped_n),
                    "comment": comment,
                })

    with open(output_csv, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=out_header)
        writer.writeheader()
        writer.writerows(rows)

    print(f"Wrote {len(rows)} rows to {output_csv}")
    print(f"  Rounds: {sorted(set(r['round'] for r in rows))}")
    print(f"  Cohorts: {sorted(set(r['cohort'] for r in rows))}")
    print(f"  Groups: {sorted(set(r['group'] for r in rows))}")
    print(f"  Conditions: {sorted(set(r['condition'] for r in rows))}")


if __name__ == "__main__":
    main()

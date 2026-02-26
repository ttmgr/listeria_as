#!/usr/bin/env python3
"""
Step 16 companion script: compile AMRFinderPlus outputs for reporting.
Purpose: create easy-to-read AMR summary tables from per-sample TSV files.

Outputs:
  amr_reads_overview.csv    — Barcode, Type, Gene Symbol, Class, Subclass, Read Count
  amr_contigs_overview.csv  — Assembler, Barcode, Type, Gene Symbol, Class, Subclass, Contig Count
  amr_gene_matrix.csv       — Barcode × Gene Symbol matrix (reads)

Usage:
    python3 16_compile_amr_overview.py <base_dir>
"""

import sys
import os
import glob
import pandas as pd
from collections import defaultdict

base_dir = sys.argv[1]
amr_dir = os.path.join(base_dir, 'processing/amrfinder')
out_dir = os.path.join(base_dir, 'processing/amrfinder/overview')
os.makedirs(out_dir, exist_ok=True)

AMR_COLS = [
    'Protein id', 'Contig id', 'Start', 'Stop', 'Strand',
    'Element symbol', 'Element name', 'Scope', 'Type', 'Subtype',
    'Class', 'Subclass', 'Method', 'Target length',
    'Reference sequence length', '% Coverage of reference',
    '% Identity to reference', 'Alignment length',
    'Closest reference accession', 'Closest reference name',
    'HMM accession', 'HMM description'
]


def parse_amr_files(subdir, source_label):
    """Parse all AMRFinder TSV files in a directory."""
    rows = []
    pattern = os.path.join(amr_dir, subdir, f'amrfinder_{subdir}_*.tsv')
    for f in sorted(glob.glob(pattern)):
        basename = os.path.basename(f)
        # Extract sample name: amrfinder_{subdir}_{SAMPLE}.tsv
        sample = basename.replace(f'amrfinder_{subdir}_', '').replace('.tsv', '')

        try:
            df = pd.read_csv(f, sep='\t', comment='#')
        except Exception:
            continue

        if df.empty or len(df) == 0:
            continue

        # Standardize column names (AMRFinder output can vary slightly)
        gene_col = 'Element symbol' if 'Element symbol' in df.columns else df.columns[5]
        class_col = 'Class' if 'Class' in df.columns else df.columns[10]
        subclass_col = 'Subclass' if 'Subclass' in df.columns else df.columns[11]

        for _, row in df.iterrows():
            rows.append({
                'source': source_label,
                'sample': sample,
                'gene': row[gene_col],
                'class': row[class_col],
                'subclass': row[subclass_col],
            })
    return rows


# ---- 1. Parse reads ----
print("Parsing AMRFinder reads results...")
reads_rows = parse_amr_files('reads', 'reads')

# ---- 2. Parse contigs (Flye + mdbg + myloasm) ----
print("Parsing AMRFinder Flye results...")
flye_rows = parse_amr_files('flye', 'flye')

print("Parsing AMRFinder mdbg results...")
mdbg_rows = parse_amr_files('mdbg', 'mdbg')

print("Parsing AMRFinder Myloasm results...")
myloasm_rows = parse_amr_files('myloasm', 'myloasm')

# ---- 3. Build reads overview ----
if reads_rows:
    df_reads = pd.DataFrame(reads_rows)
    # Extract barcode and type
    df_reads['barcode'] = df_reads['sample'].str.extract(r'(barcode\d+)')[0]
    df_reads['type'] = df_reads['sample'].str.extract(r'barcode\d+_(\w+)')[0]

    # Aggregate: barcode, type, gene → count of reads
    reads_overview = (df_reads.groupby(['barcode', 'type', 'gene', 'class', 'subclass'])
                      .size().reset_index(name='num_reads')
                      .sort_values(['barcode', 'type', 'num_reads'], ascending=[True, True, False]))

    # Rename for final output
    reads_overview = reads_overview.rename(columns={
        'barcode': 'Barcode', 
        'type': 'Type', 
        'gene': 'Gene Symbol', 
        'class': 'Class', 
        'subclass': 'Subclass', 
        'num_reads': 'Read Count'
    })

    reads_csv = os.path.join(out_dir, 'amr_reads_overview.csv')
    reads_overview.to_csv(reads_csv, index=False)
    print(f"Saved: {reads_csv} ({len(reads_overview)} rows)")

    # Gene matrix: barcode_type × gene
    matrix = df_reads.pivot_table(index='sample', columns='gene', aggfunc='size', fill_value=0)
    matrix.to_csv(os.path.join(out_dir, 'amr_gene_matrix.csv'))
else:
    print("WARNING: No reads AMRFinder results found")

# ---- 4. Build contigs overview ----
contig_rows = flye_rows + mdbg_rows + myloasm_rows
if contig_rows:
    df_contigs = pd.DataFrame(contig_rows)
    df_contigs['barcode'] = df_contigs['sample'].str.extract(r'(barcode\d+)')[0]
    df_contigs['type'] = df_contigs['sample'].str.extract(r'barcode\d+_(\w+)')[0]

    # Aggregate: assembler, barcode, type, gene → count of contigs
    contigs_overview = (df_contigs.groupby(['source', 'barcode', 'type', 'gene', 'class', 'subclass'])
                        .size().reset_index(name='num_contigs')
                        .rename(columns={'source': 'assembler'})
                        .sort_values(['assembler', 'barcode', 'type', 'num_contigs'],
                                     ascending=[True, True, True, False]))
    
    # Rename for final output
    contigs_overview = contigs_overview.rename(columns={
        'assembler': 'Assembler',
        'barcode': 'Barcode',
        'type': 'Type', 
        'gene': 'Gene Symbol', 
        'class': 'Class', 
        'subclass': 'Subclass', 
        'num_contigs': 'Contig Count'
    })

    contigs_csv = os.path.join(out_dir, 'amr_contigs_overview.csv')
    contigs_overview.to_csv(contigs_csv, index=False)
    print(f"Saved: {contigs_csv} ({len(contigs_overview)} rows)")
else:
    print("WARNING: No contigs AMRFinder results found")

# ---- 5. Print summary stats ----
print("\n" + "=" * 60)
print("  AMR OVERVIEW SUMMARY")
print("=" * 60)

if reads_rows:
    df_r = pd.DataFrame(reads_rows)
    print(f"\n  Reads:")
    print(f"    Total AMR hits:         {len(df_r):,}")
    print(f"    Unique genes:           {df_r['gene'].nunique()}")
    print(f"    Samples with AMR hits:  {df_r['sample'].nunique()}")

if contig_rows:
    df_c = pd.DataFrame(contig_rows)
    print(f"\n  Contigs:")
    print(f"    Total AMR hits:         {len(df_c):,}")
    print(f"    Unique genes:           {df_c['gene'].nunique()}")

print("\n" + "=" * 60)
print(f"All outputs saved to: {out_dir}")

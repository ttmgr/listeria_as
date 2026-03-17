#!/usr/bin/env python3
"""
Step 15 companion script: create one combined Listeria overview.
Purpose: merge read-level and contig-level Listeria results into tables and figures.

Usage:
    python3 15_compile_listeria_overview.py <base_dir>

Expects:
    <base_dir>/processing/listeria/listeria_summary.tsv
    <base_dir>/processing/listeria/listeria_contigs_summary.tsv
    <base_dir>/processing/nanostat/barcode*  (for total read counts)

Main output:
    <base_dir>/processing/listeria/overview/listeria_overview.csv
"""

import sys
import os
import re
import glob
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import numpy as np

# ---------- Config ----------
plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 300,
    'savefig.bbox': 'tight',
})

C_AS = '#2196F3'
C_N  = '#FF9800'
C_FLYE = '#4CAF50'
C_MDBG = '#9C27B0'
C_MYLOASM = '#E91E63'  # Pink/Red color for myloasm

# ---------- Paths ----------
base_dir = sys.argv[1]
listeria_tsv = os.path.join(base_dir, 'processing/listeria/listeria_summary.tsv')
contigs_tsv  = os.path.join(base_dir, 'processing/listeria/listeria_contigs_summary.tsv')
nanostat_dir = os.path.join(base_dir, 'processing/nanostat')
out_dir      = os.path.join(base_dir, 'processing/listeria/overview')
os.makedirs(out_dir, exist_ok=True)

# ---------- 1. Parse read-level Listeria data ----------
reads_cols = ['sample', 'listeria_reads', 'listeria_bases', 'mean_read_len', 'median_read_len']
df_reads = pd.read_csv(listeria_tsv, sep='\t', header=None, names=reads_cols)
for col in reads_cols[1:]:
    df_reads[col] = pd.to_numeric(df_reads[col], errors='coerce').fillna(0)

# Remove duplicates (keep last = most recent run)
df_reads = df_reads.drop_duplicates(subset='sample', keep='last')

# ---------- 2. Parse total read counts from NanoStat ----------
total_reads_data = []

# Try multiple naming patterns — NanoStat output varies
patterns = [
    os.path.join(nanostat_dir, 'nanostat_r*_barcode*'),         # round-prefixed (new naming)
    os.path.join(nanostat_dir, 'nanostat_barcode*'),            # renamed by 04_nanostat.sh (legacy)
    os.path.join(nanostat_dir, 'r*_barcode*'),                  # round-prefixed without nanostat_ prefix
    os.path.join(nanostat_dir, 'barcode*NanoStats.txt'),        # NanoStat default name
    os.path.join(nanostat_dir, 'barcode*'),                     # any barcode file (legacy)
]

stat_files = []
for pat in patterns:
    stat_files = glob.glob(pat)
    if stat_files:
        print(f"  Found {len(stat_files)} NanoStat files with pattern: {os.path.basename(pat)}")
        break

if not stat_files:
    print(f"  WARNING: No NanoStat files found in {nanostat_dir}")
    print(f"  Tried patterns: {[os.path.basename(p) for p in patterns]}")
    try:
        contents = os.listdir(nanostat_dir)[:10]
        print(f"  Directory contents (first 10): {contents}")
    except:
        print(f"  Cannot list directory")

for stat_file in stat_files:
    # Extract sample name: strip 'nanostat_' prefix and file extension
    fname = os.path.basename(stat_file)
    sample = fname.replace('nanostat_', '').replace('NanoStats.txt', '').replace('.txt', '').replace('.tsv', '')
    # Ensure we have a valid sample name (round-prefixed or legacy format)
    if not re.match(r'(r\d+_)?barcode\d+_(AS|N)', sample):
        continue
    total_reads = 0
    total_bases = 0
    with open(stat_file) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) == 2:
                if parts[0] == 'number_of_reads':
                    total_reads = int(float(parts[1]))
                elif parts[0] == 'number_of_bases':
                    total_bases = int(float(parts[1]))
    total_reads_data.append({'sample': sample, 'total_reads': total_reads, 'total_bases': total_bases})

df_total = pd.DataFrame(total_reads_data)
print(f"  NanoStat data: {len(df_total)} samples parsed")

# Fallback: if no NanoStat data, try qc_metrics from report
if len(df_total) == 0:
    qc_path = os.path.join(base_dir, 'processing/report/qc_metrics.csv')
    if os.path.exists(qc_path):
        print(f"  Falling back to qc_metrics.csv for total reads")
        df_qc = pd.read_csv(qc_path)
        df_total = df_qc.rename(columns={
            'Sample': 'sample',
            'Total Reads': 'total_reads',
            'Total Bases (bp)': 'total_bases'
        })[['sample', 'total_reads', 'total_bases']].copy()
        df_total['total_reads'] = pd.to_numeric(df_total['total_reads'], errors='coerce').fillna(0).astype(int)
        df_total['total_bases'] = pd.to_numeric(df_total['total_bases'], errors='coerce').fillna(0).astype(int)
        print(f"  Loaded {len(df_total)} samples from qc_metrics.csv")

# ---------- 3. Parse contig-level Listeria data ----------
contigs_cols = ['sample', 'assembler', 'contig_count', 'contig_bases', 'median_contig_len', 'kraken_hits']
if os.path.exists(contigs_tsv):
    df_contigs = pd.read_csv(contigs_tsv, sep='\t', header=None, names=contigs_cols)
    for col in contigs_cols[2:]:
        df_contigs[col] = pd.to_numeric(df_contigs[col], errors='coerce').fillna(0)
    df_contigs = df_contigs.drop_duplicates(subset=['sample', 'assembler'], keep='last')

    # Pivot: one row per sample, columns for each assembler
    flye = df_contigs[df_contigs['assembler'] == 'flye'][['sample', 'contig_count', 'contig_bases', 'median_contig_len']].rename(
        columns={'contig_count': 'flye_contigs', 'contig_bases': 'flye_contig_bases', 'median_contig_len': 'flye_median_contig_len'})
    mdbg = df_contigs[df_contigs['assembler'] == 'mdbg'][['sample', 'contig_count', 'contig_bases', 'median_contig_len']].rename(
        columns={'contig_count': 'mdbg_contigs', 'contig_bases': 'mdbg_contig_bases', 'median_contig_len': 'mdbg_median_contig_len'})
    myloasm = df_contigs[df_contigs['assembler'] == 'myloasm'][['sample', 'contig_count', 'contig_bases', 'median_contig_len']].rename(
        columns={'contig_count': 'myloasm_contigs', 'contig_bases': 'myloasm_contig_bases', 'median_contig_len': 'myloasm_median_contig_len'})
else:
    print("WARNING: No contigs summary found, skipping contig data")
    flye = pd.DataFrame(columns=['sample', 'flye_contigs', 'flye_contig_bases', 'flye_median_contig_len'])
    mdbg = pd.DataFrame(columns=['sample', 'mdbg_contigs', 'mdbg_contig_bases', 'mdbg_median_contig_len'])
    myloasm = pd.DataFrame(columns=['sample', 'myloasm_contigs', 'myloasm_contig_bases', 'myloasm_median_contig_len'])

# ---------- 4. Merge everything ----------
df = df_reads.merge(df_total, on='sample', how='left')
df = df.merge(flye, on='sample', how='left')
df = df.merge(mdbg, on='sample', how='left')
df = df.merge(myloasm, on='sample', how='left')
df = df.fillna(0)

# Derived columns
df['pct_listeria'] = np.where(df['total_reads'] > 0,
    (df['listeria_reads'] / df['total_reads'] * 100).round(4), 0)
df['round'] = df['sample'].str.extract(r'^(r\d+)_')[0].fillna('')
df['barcode'] = df['sample'].str.extract(r'(?:r\d+_)?(barcode\d+)')[0]
df['barcode_num'] = df['barcode'].str.extract(r'barcode(\d+)')[0].astype(int)
df['type'] = df['sample'].str.extract(r'barcode\d+_(\w+)')[0]

# Sort
df = df.sort_values(['barcode_num', 'type'])

# Ensure int columns
for col in ['listeria_reads', 'total_reads', 'flye_contigs', 'mdbg_contigs', 'myloasm_contigs']:
    df[col] = df[col].astype(int)

# ---------- 5. Save comprehensive CSV ----------
# Rename columns for final report
df_final = df.rename(columns={
    'sample': 'Sample',
    'listeria_reads': 'Listeria Reads',
    'listeria_bases': 'Listeria Bases',
    'mean_read_len': 'Mean Read Length',
    'median_read_len': 'Median Read Length',
    'total_reads': 'Total Reads',
    'total_bases': 'Total Bases',
    'flye_contigs': 'Flye Contigs',
    'mdbg_contigs': 'MetaMDBG Contigs',
    'myloasm_contigs': 'Myloasm Contigs',
    'pct_listeria': 'Listeria (%)',
    'barcode': 'Barcode',
    'type': 'Type'
})

out_csv = os.path.join(out_dir, 'listeria_overview.csv')
df_final.to_csv(out_csv, index=False)
print(f"Saved: {out_csv}")

# ---------- 6. Print formatted overview ----------
print("\n" + "=" * 80)
print("  COMPREHENSIVE LISTERIA OVERVIEW")
print("=" * 80)

for stype in ['AS', 'N']:
    sub = df[df['type'] == stype]
    print(f"\n  {'='*35} {stype} SAMPLES {'='*35}")
    print(f"  Total samples:               {len(sub)}")
    print(f"  Samples with Listeria reads: {(sub['listeria_reads'] > 0).sum()}/{len(sub)}")
    print(f"  Total Listeria reads:        {sub['listeria_reads'].sum():,}")
    print(f"  Mean Listeria reads/sample:  {sub['listeria_reads'].mean():,.1f}")
    print(f"  Median Listeria reads:       {sub['listeria_reads'].median():,.1f}")
    print(f"  Mean % Listeria:             {sub['pct_listeria'].mean():.3f}%")
    if sub['listeria_reads'].sum() > 0:
        weighted_mean = (sub['listeria_bases'].sum() / sub['listeria_reads'].sum())
        print(f"  Weighted mean read length:   {weighted_mean:,.1f} bp")
    print(f"  Flye Listeria contigs:       {sub['flye_contigs'].sum():,.0f}")
    print(f"  mdbg Listeria contigs:       {sub['mdbg_contigs'].sum():,.0f}")
    print(f"  Myloasm Listeria contigs:    {sub['myloasm_contigs'].sum():,.0f}")

# Per-sample table
print(f"\n  {'='*88}")
print(f"  {'Sample':<18} {'Reads':>8} {'%List':>7} {'MeanLen':>8} {'MedLen':>8} {'FlyeCtg':>8} {'MdbgCtg':>8} {'MyloCtg':>8}")
print(f"  {'-'*18} {'-'*8} {'-'*7} {'-'*8} {'-'*8} {'-'*8} {'-'*8} {'-'*8}")
for _, row in df.iterrows():
    print(f"  {row['sample']:<18} {int(row['listeria_reads']):>8,} {row['pct_listeria']:>6.3f}% "
          f"{row['mean_read_len']:>8.1f} {row['median_read_len']:>8.1f} "
          f"{int(row['flye_contigs']):>8} {int(row['mdbg_contigs']):>8} {int(row['myloasm_contigs']):>8}")

# ---------- 7. PLOTS ----------

barcodes = sorted(df['barcode_num'].unique())
x = np.arange(len(barcodes))
width = 0.38

# ---- Plot 1: % Listeria reads, AS vs N ----
fig, ax = plt.subplots(figsize=(16, 5))
pivot_pct = df.pivot(index='barcode_num', columns='type', values='pct_listeria').fillna(0).reindex(barcodes)
if 'AS' in pivot_pct.columns:
    ax.bar(x - width/2, pivot_pct['AS'], width, label='AS', color=C_AS, edgecolor='white', linewidth=0.5)
if 'N' in pivot_pct.columns:
    ax.bar(x + width/2, pivot_pct['N'], width, label='N', color=C_N, edgecolor='white', linewidth=0.5)
ax.set_xlabel('Barcode')
ax.set_ylabel('% Listeria reads')
ax.set_title('Percentage of reads classified as Listeria')
ax.set_xticks(x)
ax.set_xticklabels([f'{b:02d}' for b in barcodes], rotation=45, ha='right', fontsize=9)
ax.legend(title='Sample type', frameon=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
try:
    plt.tight_layout()
except Exception:
    pass
fig.savefig(os.path.join(out_dir, 'pct_listeria_per_barcode.pdf'))
fig.savefig(os.path.join(out_dir, 'pct_listeria_per_barcode.png'))
plt.close()
print(f"\nSaved: pct_listeria_per_barcode.pdf/.png")

# ---- Plot 2: Listeria read counts (log scale) ----
fig, ax = plt.subplots(figsize=(16, 5))
pivot_reads = df.pivot(index='barcode_num', columns='type', values='listeria_reads').fillna(0).reindex(barcodes)
if 'AS' in pivot_reads.columns:
    ax.bar(x - width/2, pivot_reads['AS'].replace(0, np.nan), width, label='AS', color=C_AS, edgecolor='white', linewidth=0.5)
if 'N' in pivot_reads.columns:
    ax.bar(x + width/2, pivot_reads['N'].replace(0, np.nan), width, label='N', color=C_N, edgecolor='white', linewidth=0.5)
ax.set_xlabel('Barcode')
ax.set_ylabel('Listeria read count (log)')
ax.set_title('Listeria read counts per barcode')
ax.set_xticks(x)
ax.set_xticklabels([f'{b:02d}' for b in barcodes], rotation=45, ha='right', fontsize=9)
ax.set_yscale('log')
ax.legend(title='Sample type', frameon=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
try:
    plt.tight_layout()
except Exception:
    pass
fig.savefig(os.path.join(out_dir, 'listeria_reads_log_per_barcode.pdf'))
fig.savefig(os.path.join(out_dir, 'listeria_reads_log_per_barcode.png'))
plt.close()
print(f"Saved: listeria_reads_log_per_barcode.pdf/.png")

# ---- Plot 3: Contig comparison Flye vs mdbg ----
fig, axes = plt.subplots(1, 2, figsize=(14, 5))

# 3a: Contig counts
ax = axes[0]
width3 = 0.25
pivot_flye_ctg = df.pivot(index='barcode_num', columns='type', values='flye_contigs').fillna(0).reindex(barcodes)
pivot_mdbg_ctg = df.pivot(index='barcode_num', columns='type', values='mdbg_contigs').fillna(0).reindex(barcodes)
pivot_myloasm_ctg = df.pivot(index='barcode_num', columns='type', values='myloasm_contigs').fillna(0).reindex(barcodes)
# Sum AS and N for each assembler per barcode for simplicity
flye_sum = (pivot_flye_ctg.get('AS', 0) + pivot_flye_ctg.get('N', 0)) if len(pivot_flye_ctg.columns) > 0 else pd.Series(0, index=barcodes)
mdbg_sum = (pivot_mdbg_ctg.get('AS', 0) + pivot_mdbg_ctg.get('N', 0)) if len(pivot_mdbg_ctg.columns) > 0 else pd.Series(0, index=barcodes)
myloasm_sum = (pivot_myloasm_ctg.get('AS', 0) + pivot_myloasm_ctg.get('N', 0)) if len(pivot_myloasm_ctg.columns) > 0 else pd.Series(0, index=barcodes)

ax.bar(x - width3, flye_sum, width3, label='Flye+Racon', color=C_FLYE, edgecolor='white', linewidth=0.5)
ax.bar(x, mdbg_sum, width3, label='metaMDBG', color=C_MDBG, edgecolor='white', linewidth=0.5)
ax.bar(x + width3, myloasm_sum, width3, label='Myloasm', color=C_MYLOASM, edgecolor='white', linewidth=0.5)
ax.set_xlabel('Barcode')
ax.set_ylabel('Listeria contig count')
ax.set_title('Listeria contigs: Flye vs mdbg vs Myloasm')
ax.set_xticks(x)
ax.set_xticklabels([f'{b:02d}' for b in barcodes], rotation=45, ha='right', fontsize=8)
ax.legend(frameon=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# 3b: Median contig length
ax = axes[1]
pivot_flye_med = df.pivot(index='barcode_num', columns='type', values='flye_median_contig_len').fillna(0).reindex(barcodes)
pivot_mdbg_med = df.pivot(index='barcode_num', columns='type', values='mdbg_median_contig_len').fillna(0).reindex(barcodes)
pivot_myloasm_med = df.pivot(index='barcode_num', columns='type', values='myloasm_median_contig_len').fillna(0).reindex(barcodes)
flye_med = pivot_flye_med.max(axis=1).replace(0, np.nan)
mdbg_med = pivot_mdbg_med.max(axis=1).replace(0, np.nan)
myloasm_med = pivot_myloasm_med.max(axis=1).replace(0, np.nan)
ax.bar(x - width3, flye_med, width3, label='Flye+Racon', color=C_FLYE, edgecolor='white', linewidth=0.5)
ax.bar(x, mdbg_med, width3, label='metaMDBG', color=C_MDBG, edgecolor='white', linewidth=0.5)
ax.bar(x + width3, myloasm_med, width3, label='Myloasm', color=C_MYLOASM, edgecolor='white', linewidth=0.5)
ax.set_xlabel('Barcode')
ax.set_ylabel('Median Listeria contig length (bp)')
ax.set_title('Median contig length: Flye vs mdbg vs Myloasm')
ax.set_xticks(x)
ax.set_xticklabels([f'{b:02d}' for b in barcodes], rotation=45, ha='right', fontsize=8)
ax.legend(frameon=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

try:
    plt.tight_layout()
except Exception:
    pass
fig.savefig(os.path.join(out_dir, 'listeria_contigs_comparison.pdf'))
fig.savefig(os.path.join(out_dir, 'listeria_contigs_comparison.png'))
plt.close()
print(f"Saved: listeria_contigs_comparison.pdf/.png")

# ---- Plot 4: Multi-panel overview (reads + contigs + %) ----
fig, axes = plt.subplots(3, 1, figsize=(16, 12), sharex=True)

# Panel 1: Read counts
ax = axes[0]
if 'AS' in pivot_reads.columns:
    ax.bar(x - width/2, pivot_reads['AS'].replace(0, np.nan), width, label='AS', color=C_AS, edgecolor='white', linewidth=0.5)
if 'N' in pivot_reads.columns:
    ax.bar(x + width/2, pivot_reads['N'].replace(0, np.nan), width, label='N', color=C_N, edgecolor='white', linewidth=0.5)
ax.set_ylabel('Listeria reads')
ax.set_title('Listeria reads, contigs, and enrichment across barcodes')
ax.set_yscale('symlog', linthresh=1)
ax.legend(title='Type', frameon=False, loc='upper right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel 2: Contig counts (combined AS + N)
ax = axes[1]
ax.bar(x - width3, flye_sum, width3, label='Flye+Racon', color=C_FLYE, edgecolor='white', linewidth=0.5)
ax.bar(x, mdbg_sum, width3, label='metaMDBG', color=C_MDBG, edgecolor='white', linewidth=0.5)
ax.bar(x + width3, myloasm_sum, width3, label='Myloasm', color=C_MYLOASM, edgecolor='white', linewidth=0.5)
ax.set_ylabel('Listeria contigs')
ax.legend(title='Assembler', frameon=False, loc='upper right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

# Panel 3: % Listeria
ax = axes[2]
if 'AS' in pivot_pct.columns:
    ax.bar(x - width/2, pivot_pct['AS'], width, label='AS', color=C_AS, edgecolor='white', linewidth=0.5)
if 'N' in pivot_pct.columns:
    ax.bar(x + width/2, pivot_pct['N'], width, label='N', color=C_N, edgecolor='white', linewidth=0.5)
ax.set_xlabel('Barcode')
ax.set_ylabel('% Listeria')
ax.legend(title='Type', frameon=False, loc='upper right')
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

ax.set_xticks(x)
ax.set_xticklabels([f'{b:02d}' for b in barcodes], rotation=45, ha='right', fontsize=9)

try:
    plt.tight_layout()
except Exception:
    pass
fig.savefig(os.path.join(out_dir, 'listeria_multi_panel.pdf'))
fig.savefig(os.path.join(out_dir, 'listeria_multi_panel.png'))
plt.close()
print(f"Saved: listeria_multi_panel.pdf/.png")

print(f"\nAll outputs saved to: {out_dir}")
print("=" * 80)

#!/usr/bin/env python3
"""
Listeria plotting helper.
Purpose: make quick visual checks of Listeria read enrichment by sample.

Reads listeria_summary.tsv and produces:
  1. Grouped bar chart: AS vs N per barcode (read counts)
  2. Grouped bar chart: AS vs N per barcode (total bases)
  3. Paired dot plot: AS vs N comparison
  4. Summary statistics printed to stdout

Usage:
    python3 plot_listeria.py <listeria_summary.tsv> <output_dir>
"""

import sys
import os
import pandas as pd
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

COLORS = {'AS': '#2196F3', 'N': '#FF9800'}

# ---------- Parse ----------
tsv_path = sys.argv[1]
out_dir = sys.argv[2] if len(sys.argv) > 2 else '.'
os.makedirs(out_dir, exist_ok=True)

df = pd.read_csv(tsv_path, sep='\t', header=None,
                 names=['sample', 'listeria_reads', 'listeria_bases', 'listeria_mean_len'])

# Fill NaN with 0
df = df.fillna(0)
for col in ['listeria_reads', 'listeria_bases', 'listeria_mean_len']:
    df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0).astype(int)

# Extract barcode number and sample type (AS or N)
df['barcode'] = df['sample'].str.extract(r'(barcode\d+)')[0]
df['barcode_num'] = df['barcode'].str.extract(r'barcode(\d+)')[0].astype(int)
df['type'] = df['sample'].str.extract(r'barcode\d+_(\w+)')[0]

# Sort by barcode number
df = df.sort_values('barcode_num')

# Pivot for grouped bar charts
pivot_reads = df.pivot(index='barcode_num', columns='type', values='listeria_reads').fillna(0)
pivot_bases = df.pivot(index='barcode_num', columns='type', values='listeria_bases').fillna(0)

barcodes = pivot_reads.index.values
x = np.arange(len(barcodes))
width = 0.38

# ============================================================
# PLOT 1: Grouped bar chart — Listeria read counts per barcode
# ============================================================
fig, ax = plt.subplots(figsize=(16, 5))

if 'AS' in pivot_reads.columns:
    ax.bar(x - width/2, pivot_reads['AS'], width, label='AS', color=COLORS['AS'], edgecolor='white', linewidth=0.5)
if 'N' in pivot_reads.columns:
    ax.bar(x + width/2, pivot_reads['N'], width, label='N', color=COLORS['N'], edgecolor='white', linewidth=0.5)

ax.set_xlabel('Barcode')
ax.set_ylabel('Listeria read count')
ax.set_title('Listeria reads per barcode (AS vs N)')
ax.set_xticks(x)
ax.set_xticklabels([f'{b:02d}' for b in barcodes], rotation=45, ha='right', fontsize=9)
ax.legend(title='Sample type', frameon=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.yaxis.set_major_formatter(ticker.FuncFormatter(lambda v, _: f'{int(v):,}'))

plt.tight_layout()
fig.savefig(os.path.join(out_dir, 'listeria_reads_per_barcode.pdf'))
fig.savefig(os.path.join(out_dir, 'listeria_reads_per_barcode.png'))
plt.close()
print(f"Saved: listeria_reads_per_barcode.pdf/.png")

# ============================================================
# PLOT 2: Grouped bar chart — Listeria total bases per barcode
# ============================================================
fig, ax = plt.subplots(figsize=(16, 5))

if 'AS' in pivot_bases.columns:
    ax.bar(x - width/2, pivot_bases['AS'] / 1e6, width, label='AS', color=COLORS['AS'], edgecolor='white', linewidth=0.5)
if 'N' in pivot_bases.columns:
    ax.bar(x + width/2, pivot_bases['N'] / 1e6, width, label='N', color=COLORS['N'], edgecolor='white', linewidth=0.5)

ax.set_xlabel('Barcode')
ax.set_ylabel('Listeria bases (Mb)')
ax.set_title('Listeria total bases per barcode (AS vs N)')
ax.set_xticks(x)
ax.set_xticklabels([f'{b:02d}' for b in barcodes], rotation=45, ha='right', fontsize=9)
ax.legend(title='Sample type', frameon=False)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(os.path.join(out_dir, 'listeria_bases_per_barcode.pdf'))
fig.savefig(os.path.join(out_dir, 'listeria_bases_per_barcode.png'))
plt.close()
print(f"Saved: listeria_bases_per_barcode.pdf/.png")

# ============================================================
# PLOT 3: Paired dot plot — AS vs N per barcode
# ============================================================
fig, ax = plt.subplots(figsize=(8, 6))

for bc in barcodes:
    row_as = df[(df['barcode_num'] == bc) & (df['type'] == 'AS')]
    row_n = df[(df['barcode_num'] == bc) & (df['type'] == 'N')]
    val_as = row_as['listeria_reads'].values[0] if len(row_as) > 0 else 0
    val_n = row_n['listeria_reads'].values[0] if len(row_n) > 0 else 0

    # Only draw line if at least one is non-zero
    if val_as > 0 or val_n > 0:
        ax.plot([0, 1], [val_as, val_n], '-', color='#999999', alpha=0.4, linewidth=0.8)

# Plot all AS and N points
as_vals = []
n_vals = []
for bc in barcodes:
    row_as = df[(df['barcode_num'] == bc) & (df['type'] == 'AS')]
    row_n = df[(df['barcode_num'] == bc) & (df['type'] == 'N')]
    as_vals.append(row_as['listeria_reads'].values[0] if len(row_as) > 0 else 0)
    n_vals.append(row_n['listeria_reads'].values[0] if len(row_n) > 0 else 0)

ax.scatter([0]*len(as_vals), as_vals, color=COLORS['AS'], s=40, zorder=5, label='AS', alpha=0.8)
ax.scatter([1]*len(n_vals), n_vals, color=COLORS['N'], s=40, zorder=5, label='N', alpha=0.8)

ax.set_xticks([0, 1])
ax.set_xticklabels(['AS', 'N'], fontsize=13)
ax.set_ylabel('Listeria read count')
ax.set_title('Listeria reads: AS vs N (paired by barcode)')
ax.set_yscale('symlog', linthresh=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
ax.legend(frameon=False)

plt.tight_layout()
fig.savefig(os.path.join(out_dir, 'listeria_AS_vs_N_paired.pdf'))
fig.savefig(os.path.join(out_dir, 'listeria_AS_vs_N_paired.png'))
plt.close()
print(f"Saved: listeria_AS_vs_N_paired.pdf/.png")

# ============================================================
# PLOT 4: Box/strip plot — AS vs N distribution
# ============================================================
fig, ax = plt.subplots(figsize=(5, 6))

as_data = df[df['type'] == 'AS']['listeria_reads'].values
n_data = df[df['type'] == 'N']['listeria_reads'].values

bp = ax.boxplot([as_data, n_data], positions=[0, 1], widths=0.5,
                patch_artist=True, showfliers=False,
                medianprops=dict(color='black', linewidth=1.5))

bp['boxes'][0].set_facecolor(COLORS['AS'])
bp['boxes'][0].set_alpha(0.3)
bp['boxes'][1].set_facecolor(COLORS['N'])
bp['boxes'][1].set_alpha(0.3)

# Overlay strip plot with jitter
np.random.seed(42)
jitter_as = np.random.normal(0, 0.04, len(as_data))
jitter_n = np.random.normal(1, 0.04, len(n_data))
ax.scatter(jitter_as, as_data, color=COLORS['AS'], s=30, alpha=0.7, zorder=5, edgecolors='white', linewidth=0.5)
ax.scatter(jitter_n, n_data, color=COLORS['N'], s=30, alpha=0.7, zorder=5, edgecolors='white', linewidth=0.5)

ax.set_xticks([0, 1])
ax.set_xticklabels(['AS', 'N'], fontsize=13)
ax.set_ylabel('Listeria read count')
ax.set_title('Listeria reads: AS vs N')
ax.set_yscale('symlog', linthresh=1)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)

plt.tight_layout()
fig.savefig(os.path.join(out_dir, 'listeria_AS_vs_N_boxplot.pdf'))
fig.savefig(os.path.join(out_dir, 'listeria_AS_vs_N_boxplot.png'))
plt.close()
print(f"Saved: listeria_AS_vs_N_boxplot.pdf/.png")

# ============================================================
# Summary statistics
# ============================================================
print("\n" + "="*50)
print("  LISTERIA SUMMARY STATISTICS")
print("="*50)
for stype in ['AS', 'N']:
    subset = df[df['type'] == stype]
    reads = subset['listeria_reads']
    print(f"\n  {stype} samples (n={len(subset)}):")
    print(f"    Total Listeria reads:  {reads.sum():,}")
    print(f"    Mean per sample:       {reads.mean():,.1f}")
    print(f"    Median per sample:     {reads.median():,.1f}")
    print(f"    Max:                   {reads.max():,}")
    print(f"    Samples with >0:       {(reads > 0).sum()}/{len(subset)}")

# Fold change
as_total = df[df['type'] == 'AS']['listeria_reads'].sum()
n_total = df[df['type'] == 'N']['listeria_reads'].sum()
if n_total > 0:
    print(f"\n  AS/N fold enrichment:    {as_total/n_total:.1f}x")
else:
    print(f"\n  AS/N fold enrichment:    ∞ (no N reads)")

print(f"\n  Output directory: {out_dir}")
print("="*50)

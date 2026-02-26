#!/usr/bin/env python3
"""
Step 20 comparison report.
Purpose: compare Adaptive Sampling (AS) vs Normal (N) across the defined Black group.

Usage:
    python3 20_comparison_report.py <base_dir>

Expects:
    <base_dir>/sample_metadata.csv
    <base_dir>/processing/stats/read_metrics_summary.csv
    <base_dir>/processing/listeria/overview/listeria_overview.csv
    <base_dir>/processing/amrfinder/overview/amr_reads_overview.csv
    <base_dir>/processing/amrfinder/overview/amr_contigs_overview.csv
    <base_dir>/processing/stats/assembly_stats_{flye,mdbg,myloasm}.tsv
"""

import sys
import os
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io
import base64

# ============================================================
# Config
# ============================================================

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.sans-serif': ['Arial', 'DejaVu Sans'],
    'font.size': 11,
    'axes.titlesize': 13,
    'axes.labelsize': 12,
    'figure.dpi': 150,
    'savefig.dpi': 200,
    'savefig.bbox': 'tight',
})

# Colors
C_AS = '#86efac'      # Green (Adaptive Sampling)
C_N = '#fca5a5'       # Red (Normal)
C_GROUP_A = '#3b82f6'  # Blue  — Sponge / PowerSoil
C_GROUP_C = '#f59e0b'  # Amber — Cotton / Zymo
C_GROUP_D = '#8b5cf6'  # Purple — Zymo / Zymo
GROUP_COLORS = {'Black_A': C_GROUP_A, 'Black_C': C_GROUP_C, 'Black_D': C_GROUP_D}

# Black barcodes only
BLACK_BARCODES = [3, 4, 5, 6, 7, 14, 15, 16, 18, 19, 26, 27, 28, 29, 30]

base_dir = sys.argv[1]
out_dir = os.path.join(base_dir, 'processing/report')
os.makedirs(out_dir, exist_ok=True)

# ============================================================
# Helpers
# ============================================================

def fig_to_b64(fig):
    """Convert matplotlib figure to base64 PNG string."""
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')


def df_to_html(df, table_id=''):
    """Convert DataFrame to styled HTML table."""
    html = f'<table id="{table_id}">\n<thead><tr>'
    for col in df.columns:
        html += f'<th onclick="sortTable(this)">{col}</th>'
    html += '</tr></thead>\n<tbody>\n'
    for _, row in df.iterrows():
        html += '<tr>'
        for val in row:
            if isinstance(val, float):
                html += f'<td>{val:,.2f}</td>'
            elif isinstance(val, (int, np.integer)):
                html += f'<td>{val:,}</td>'
            else:
                html += f'<td>{val}</td>'
        html += '</tr>\n'
    html += '</tbody>\n</table>'
    return html


# ============================================================
# 1. Load Metadata
# ============================================================

print("Loading sample metadata...")
meta_path = os.path.join(base_dir, 'sample_metadata.csv')
df_meta = pd.read_csv(meta_path)

# Filter to Black samples only
df_black = df_meta[df_meta['barcode'].isin(BLACK_BARCODES)].copy()
# Create filename key: barcode03_AS, barcode03_N, etc.
df_black['sample'] = df_black.apply(
    lambda r: f"barcode{r['barcode']:02d}_{r['condition']}", axis=1)

print(f"  Black samples: {len(df_black)} rows ({len(df_black)//2} barcodes × 2 conditions)")

# ============================================================
# 2. Load Pipeline Outputs
# ============================================================

print("Loading pipeline outputs...")

# QC metrics
qc_path = os.path.join(base_dir, 'processing/stats/read_metrics_summary.csv')
df_qc = pd.read_csv(qc_path) if os.path.exists(qc_path) else pd.DataFrame()

# Listeria overview
listeria_path = os.path.join(base_dir, 'processing/listeria/overview/listeria_overview.csv')
df_listeria = pd.read_csv(listeria_path) if os.path.exists(listeria_path) else pd.DataFrame()

# AMR
amr_reads_path = os.path.join(base_dir, 'processing/amrfinder/overview/amr_reads_overview.csv')
amr_contigs_path = os.path.join(base_dir, 'processing/amrfinder/overview/amr_contigs_overview.csv')
df_amr_reads = pd.read_csv(amr_reads_path) if os.path.exists(amr_reads_path) else pd.DataFrame()
df_amr_contigs = pd.read_csv(amr_contigs_path) if os.path.exists(amr_contigs_path) else pd.DataFrame()

# Assembly stats
def load_asm_stats(name):
    p = os.path.join(base_dir, f'processing/stats/assembly_stats_{name}.tsv')
    if os.path.exists(p):
        return pd.read_csv(p, sep='\t')
    return pd.DataFrame()

df_asm_flye = load_asm_stats('flye')
df_asm_mdbg = load_asm_stats('mdbg')
df_asm_myloasm = load_asm_stats('myloasm')

# ============================================================
# 3. Merge Metadata with QC
# ============================================================

print("Merging data...")

# Standardize sample column names across datasets
def standardize_sample_col(df):
    """Find the sample column and rename to 'sample'."""
    for col in ['Sample', 'sample', 'Barcode']:
        if col in df.columns and col != 'sample':
            df = df.rename(columns={col: 'sample'})
            break
    return df

df_qc = standardize_sample_col(df_qc)
df_listeria = standardize_sample_col(df_listeria)

# Merge QC with metadata
if len(df_qc) > 0 and 'sample' in df_qc.columns:
    df = df_black.merge(df_qc, on='sample', how='left')
else:
    df = df_black.copy()

# Merge Listeria data
if len(df_listeria) > 0 and 'sample' in df_listeria.columns:
    # Avoid duplicate columns
    listeria_cols = [c for c in df_listeria.columns if c not in df.columns or c == 'sample']
    df = df.merge(df_listeria[listeria_cols], on='sample', how='left')

df = df.fillna(0)

# Ensure numeric columns
for col in ['dna_concentration_ng_ul']:
    if col in df.columns:
        df[col] = pd.to_numeric(df[col], errors='coerce').fillna(0)

# ============================================================
# 4. Compute Enrichment Ratios
# ============================================================

print("Computing enrichment ratios...")

# Get listeria reads column
listeria_col = None
for candidate in ['listeria_reads', 'Listeria Reads']:
    if candidate in df.columns:
        listeria_col = candidate
        break

if listeria_col:
    df[listeria_col] = pd.to_numeric(df[listeria_col], errors='coerce').fillna(0)

# Compute AS/N enrichment per barcode
enrichment_rows = []
for bc in BLACK_BARCODES:
    as_row = df[(df['barcode'] == bc) & (df['condition'] == 'AS')]
    n_row = df[(df['barcode'] == bc) & (df['condition'] == 'N')]
    if len(as_row) > 0 and len(n_row) > 0:
        as_reads = as_row[listeria_col].values[0] if listeria_col else 0
        n_reads = n_row[listeria_col].values[0] if listeria_col else 0
        abs_diff = as_reads - n_reads
        ratio = (as_reads / n_reads) if n_reads > 0 else float('inf') if as_reads > 0 else 1.0

        # Get total reads for relative enrichment
        as_total = 0
        n_total = 0
        for rc in ['number_of_reads', 'Total Reads', 'total_reads']:
            if rc in as_row.columns:
                as_total = float(as_row[rc].values[0])
                n_total = float(n_row[rc].values[0])
                break
        as_pct = (as_reads / as_total * 100) if as_total > 0 else 0
        n_pct = (n_reads / n_total * 100) if n_total > 0 else 0
        relative_enrichment = (as_pct / n_pct) if n_pct > 0 else float('inf') if as_pct > 0 else 1.0

        # Get read length stats
        as_mean_len = 0
        n_mean_len = 0
        as_median_len = 0
        n_median_len = 0
        for ml in ['listeria_mean_len', 'Mean Read Length', 'mean_read_len']:
            if ml in as_row.columns:
                as_mean_len = float(as_row[ml].values[0])
                n_mean_len = float(n_row[ml].values[0])
                break
        for mdl in ['listeria_median_len', 'Median Read Length', 'median_read_len']:
            if mdl in as_row.columns:
                as_median_len = float(as_row[mdl].values[0])
                n_median_len = float(n_row[mdl].values[0])
                break
        enrichment_rows.append({
            'barcode': bc,
            'sample_id': as_row['sample_id'].values[0],
            'group': as_row['group'].values[0],
            'swab_type': as_row['swab_type'].values[0],
            'kit': as_row['kit'].values[0],
            'dna_ng_ul': as_row['dna_concentration_ng_ul'].values[0],
            'listeria_AS': int(as_reads),
            'listeria_N': int(n_reads),
            'abs_difference': int(abs_diff),
            'pct_AS': round(as_pct, 4),
            'pct_N': round(n_pct, 4),
            'absolute_enrichment': round(ratio, 2),
            'relative_enrichment': round(relative_enrichment, 2),
            'as_total_reads': int(as_total),
            'n_total_reads': int(n_total),
            'as_mean_len': round(as_mean_len, 1),
            'n_mean_len': round(n_mean_len, 1),
            'as_median_len': round(as_median_len, 1),
            'n_median_len': round(n_median_len, 1),
        })

df_enrichment = pd.DataFrame(enrichment_rows)

# ============================================================
# 5. Export merged CSV
# ============================================================

export_csv = os.path.join(out_dir, 'comparison_data.csv')
df.to_csv(export_csv, index=False)
print(f"Exported: {export_csv}")

enrichment_csv = os.path.join(out_dir, 'enrichment_ratios.csv')
if len(df_enrichment) > 0:
    df_enrichment.to_csv(enrichment_csv, index=False)
    print(f"Exported: {enrichment_csv}")

# ============================================================
# 6. Generate Plots
# ============================================================

print("Generating plots...")
plots = {}

# --- Plot 1: Total reads AS vs N, paired by sample ---
fig, ax = plt.subplots(figsize=(14, 5))
reads_col = None
for candidate in ['number_of_reads', 'Total Reads', 'total_reads']:
    if candidate in df.columns:
        reads_col = candidate
        break

if reads_col:
    df_sorted = df.sort_values(['group', 'barcode', 'condition'])
    barcodes_sorted = sorted(df['barcode'].unique())
    x = np.arange(len(barcodes_sorted))
    width = 0.35

    as_vals = []
    n_vals = []
    labels = []
    for bc in barcodes_sorted:
        sub = df[df['barcode'] == bc]
        labels.append(sub['sample_id'].values[0])
        as_v = sub[sub['condition'] == 'AS'][reads_col].values
        n_v = sub[sub['condition'] == 'N'][reads_col].values
        as_vals.append(int(as_v[0]) if len(as_v) > 0 else 0)
        n_vals.append(int(n_v[0]) if len(n_v) > 0 else 0)

    ax.bar(x - width/2, as_vals, width, label='Adaptive Sampling', color=C_AS, edgecolor='white')
    ax.bar(x + width/2, n_vals, width, label='Normal', color=C_N, edgecolor='white')
    ax.set_xlabel('Sample')
    ax.set_ylabel('Total Reads')
    ax.set_title('Total Read Counts: Adaptive Sampling vs Normal')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax.legend(frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Add group shading
    group_bounds = {}
    for i, bc in enumerate(barcodes_sorted):
        grp = df[df['barcode'] == bc]['group'].values[0]
        if grp not in group_bounds:
            group_bounds[grp] = [i, i]
        group_bounds[grp][1] = i
    for grp, (lo, hi) in group_bounds.items():
        ax.axvspan(lo - 0.5, hi + 0.5, alpha=0.08, color=GROUP_COLORS.get(grp, '#ccc'))
        ax.text((lo + hi) / 2, ax.get_ylim()[1] * 0.95, grp.replace('Black_', ''),
                ha='center', fontsize=10, fontweight='bold', color=GROUP_COLORS.get(grp, '#666'))

    plt.tight_layout()
    plots['total_reads'] = fig_to_b64(fig)

# --- Plot 2: % Listeria reads, AS vs N ---
pct_col = None
for candidate in ['listeria_ratio', 'Listeria (%)', 'pct_listeria']:
    if candidate in df.columns:
        pct_col = candidate
        break

if pct_col:
    fig, ax = plt.subplots(figsize=(14, 5))
    barcodes_sorted = sorted(df['barcode'].unique())
    x = np.arange(len(barcodes_sorted))

    as_pct = []
    n_pct = []
    labels = []
    for bc in barcodes_sorted:
        sub = df[df['barcode'] == bc]
        labels.append(sub['sample_id'].values[0])
        as_v = sub[sub['condition'] == 'AS'][pct_col].values
        n_v = sub[sub['condition'] == 'N'][pct_col].values
        as_pct.append(float(as_v[0]) if len(as_v) > 0 else 0)
        n_pct.append(float(n_v[0]) if len(n_v) > 0 else 0)

    ax.bar(x - width/2, as_pct, width, label='Adaptive Sampling', color=C_AS, edgecolor='white')
    ax.bar(x + width/2, n_pct, width, label='Normal', color=C_N, edgecolor='white')
    ax.set_xlabel('Sample')
    ax.set_ylabel('% Listeria Reads')
    ax.set_title('Listeria Read Fraction: Adaptive Sampling vs Normal')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax.legend(frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plots['pct_listeria'] = fig_to_b64(fig)

# --- Plot 3: Enrichment ratio per sample ---
if len(df_enrichment) > 0:
    fig, ax = plt.subplots(figsize=(12, 5))
    df_e = df_enrichment.sort_values(['group', 'barcode'])
    bar_colors = [GROUP_COLORS.get(g, '#999') for g in df_e['group']]
    bars = ax.bar(range(len(df_e)), df_e['enrichment_ratio'], color=bar_colors, edgecolor='white')
    ax.axhline(y=1.0, color='#94a3b8', linestyle='--', linewidth=1, label='No enrichment (1×)')
    ax.set_xlabel('Sample')
    ax.set_ylabel('Enrichment Ratio (AS / N)')
    ax.set_title('Adaptive Sampling Enrichment Factor for Listeria')
    ax.set_xticks(range(len(df_e)))
    ax.set_xticklabels(df_e['sample_id'], rotation=45, ha='right', fontsize=9)
    ax.legend(frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)

    # Cap infinite values for display
    for i, (_, row) in enumerate(df_e.iterrows()):
        if row['enrichment_ratio'] > 100:
            ax.text(i, min(row['enrichment_ratio'], ax.get_ylim()[1] * 0.9),
                    f"{row['enrichment_ratio']:.0f}×", ha='center', va='bottom', fontsize=8)

    plt.tight_layout()
    plots['enrichment'] = fig_to_b64(fig)

# --- Plot 3b: Absolute vs Relative Enrichment side-by-side ---
if len(df_enrichment) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(14, 5))
    df_e = df_enrichment.sort_values(['group', 'barcode'])

    # Left: Absolute enrichment (AS reads / N reads)
    bar_colors = [GROUP_COLORS.get(g, '#999') for g in df_e['group']]
    axes[0].bar(range(len(df_e)), df_e['absolute_enrichment'], color=bar_colors, edgecolor='white')
    axes[0].axhline(y=1.0, color='#94a3b8', linestyle='--', linewidth=1)
    axes[0].set_ylabel('Absolute Enrichment (AS / N reads)')
    axes[0].set_title('Absolute Enrichment')
    axes[0].set_xticks(range(len(df_e)))
    axes[0].set_xticklabels(df_e['sample_id'], rotation=45, ha='right', fontsize=8)
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)

    # Right: Relative enrichment (% Listeria AS / % Listeria N)
    axes[1].bar(range(len(df_e)), df_e['relative_enrichment'], color=bar_colors, edgecolor='white')
    axes[1].axhline(y=1.0, color='#94a3b8', linestyle='--', linewidth=1)
    axes[1].set_ylabel('Relative Enrichment (% AS / % N)')
    axes[1].set_title('Relative Enrichment')
    axes[1].set_xticks(range(len(df_e)))
    axes[1].set_xticklabels(df_e['sample_id'], rotation=45, ha='right', fontsize=8)
    axes[1].spines['top'].set_visible(False)
    axes[1].spines['right'].set_visible(False)

    plt.tight_layout()
    plots['abs_vs_rel'] = fig_to_b64(fig)

# --- Plot 3c: Absolute Listeria read counts AS vs N ---
if listeria_col:
    fig, ax = plt.subplots(figsize=(14, 5))
    barcodes_sorted = sorted(df['barcode'].unique())
    x = np.arange(len(barcodes_sorted))
    width = 0.35
    as_lis = []
    n_lis = []
    labels = []
    for bc in barcodes_sorted:
        sub = df[df['barcode'] == bc]
        labels.append(sub['sample_id'].values[0])
        as_v = sub[sub['condition'] == 'AS'][listeria_col].values
        n_v = sub[sub['condition'] == 'N'][listeria_col].values
        as_lis.append(int(as_v[0]) if len(as_v) > 0 else 0)
        n_lis.append(int(n_v[0]) if len(n_v) > 0 else 0)

    ax.bar(x - width/2, as_lis, width, label='Adaptive Sampling', color=C_AS, edgecolor='white')
    ax.bar(x + width/2, n_lis, width, label='Normal', color=C_N, edgecolor='white')
    ax.set_xlabel('Sample')
    ax.set_ylabel('Listeria Read Count')
    ax.set_title('Absolute Listeria Read Counts: AS vs Normal')
    ax.set_xticks(x)
    ax.set_xticklabels(labels, rotation=45, ha='right', fontsize=9)
    ax.legend(frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plots['abs_listeria'] = fig_to_b64(fig)

# --- Plot 3d: Mean enrichment per group (bar chart) ---
if len(df_enrichment) > 0:
    fig, axes = plt.subplots(1, 2, figsize=(10, 4))
    for ax_i, (metric, title) in enumerate([('absolute_enrichment', 'Absolute'), ('relative_enrichment', 'Relative')]):
        means = df_enrichment.groupby('group')[metric].apply(lambda x: x.replace([np.inf], np.nan).mean())
        sds = df_enrichment.groupby('group')[metric].apply(lambda x: x.replace([np.inf], np.nan).std())
        grps = ['Black_A', 'Black_C', 'Black_D']
        vals = [means.get(g, 0) for g in grps]
        errs = [sds.get(g, 0) for g in grps]
        colors = [GROUP_COLORS[g] for g in grps]
        bars = axes[ax_i].bar(range(3), vals, yerr=errs, color=colors, edgecolor='white',
                               capsize=5, error_kw={'linewidth': 1.5})
        axes[ax_i].axhline(y=1.0, color='#94a3b8', linestyle='--', linewidth=1)
        axes[ax_i].set_xticks(range(3))
        axes[ax_i].set_xticklabels(['A (Sponge/PS)', 'C (Cotton/Zymo)', 'D (Zymo/Zymo)'], fontsize=8)
        axes[ax_i].set_ylabel(f'{title} Enrichment (×)')
        axes[ax_i].set_title(f'Mean {title} Enrichment by Group')
        axes[ax_i].spines['top'].set_visible(False)
        axes[ax_i].spines['right'].set_visible(False)
        for i, v in enumerate(vals):
            if not np.isnan(v):
                axes[ax_i].text(i, v + errs[i] * 0.1, f'{v:.1f}×', ha='center', va='bottom', fontsize=9, fontweight='bold')
    plt.tight_layout()
    plots['group_enrichment'] = fig_to_b64(fig)

# --- Plot 4: DNA concentration vs total reads ---
if reads_col and 'dna_concentration_ng_ul' in df.columns:
    fig, ax = plt.subplots(figsize=(8, 6))
    for grp, color in GROUP_COLORS.items():
        sub = df[df['group'] == grp]
        ax.scatter(sub['dna_concentration_ng_ul'], sub[reads_col],
                   c=color, label=grp.replace('Black_', 'Group '),
                   s=60, alpha=0.8, edgecolors='white', linewidths=0.5)
        # Mark AS vs N
        for _, row in sub.iterrows():
            marker = '▲' if row['condition'] == 'AS' else '●'
            ax.annotate(marker, (row['dna_concentration_ng_ul'], row[reads_col]),
                       fontsize=6, ha='center', va='center')

    ax.set_xlabel('DNA Concentration (ng/µL)')
    ax.set_ylabel('Total Reads')
    ax.set_title('DNA Yield vs Sequencing Output')
    ax.legend(frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plots['dna_vs_reads'] = fig_to_b64(fig)

# --- Plot 5: Extraction method comparison (boxplot) ---
if reads_col:
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))

    # Total reads by group
    for grp, color in GROUP_COLORS.items():
        sub = df[df['group'] == grp]
        label = grp.replace('Black_', '')
        axes[0].boxplot([sub[reads_col].values], positions=[list(GROUP_COLORS.keys()).index(grp)],
                       widths=0.6, patch_artist=True,
                       boxprops=dict(facecolor=color, alpha=0.7),
                       medianprops=dict(color='black', linewidth=2))
    axes[0].set_xticks(range(len(GROUP_COLORS)))
    axes[0].set_xticklabels([g.replace('Black_', '') + f'\n({df_black[df_black["group"]==g]["swab_type"].values[0]}/{df_black[df_black["group"]==g]["kit"].values[0]})'
                             for g in GROUP_COLORS.keys()], fontsize=9)
    axes[0].set_ylabel('Total Reads')
    axes[0].set_title('Read Yield by Extraction Method')
    axes[0].spines['top'].set_visible(False)
    axes[0].spines['right'].set_visible(False)

    # Listeria reads by group
    if listeria_col:
        for grp, color in GROUP_COLORS.items():
            sub = df[df['group'] == grp]
            axes[1].boxplot([sub[listeria_col].values], positions=[list(GROUP_COLORS.keys()).index(grp)],
                           widths=0.6, patch_artist=True,
                           boxprops=dict(facecolor=color, alpha=0.7),
                           medianprops=dict(color='black', linewidth=2))
        axes[1].set_xticks(range(len(GROUP_COLORS)))
        axes[1].set_xticklabels([g.replace('Black_', '') + f'\n({df_black[df_black["group"]==g]["swab_type"].values[0]}/{df_black[df_black["group"]==g]["kit"].values[0]})'
                                 for g in GROUP_COLORS.keys()], fontsize=9)
        axes[1].set_ylabel('Listeria Reads')
        axes[1].set_title('Listeria Detection by Extraction Method')
        axes[1].spines['top'].set_visible(False)
        axes[1].spines['right'].set_visible(False)

    plt.tight_layout()
    plots['extraction_comparison'] = fig_to_b64(fig)


# ============================================================
# 7. Build HTML Report
# ============================================================

print("Building HTML report...")

css = """
body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 0; background: #f8fafc; color: #334155; }
.container { max-width: 1200px; margin: 0 auto; padding: 20px; background: #fff; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
h1, h2, h3 { color: #1e293b; border-bottom: 2px solid #e2e8f0; padding-bottom: 10px; margin-top: 2rem; }
.header { text-align: center; padding: 2rem; background: linear-gradient(135deg, #1e293b 0%, #334155 100%); color: white; border-radius: 8px; margin-bottom: 2rem; }
.header h1 { color: white; border: none; margin: 0; }
.header p { color: #94a3b8; margin: 0.5rem 0 0; }
.summary-stats { display: flex; gap: 15px; flex-wrap: wrap; justify-content: center; margin: 1.5rem 0; }
.stat-card { flex: 1; min-width: 130px; max-width: 180px; background: #fff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 1rem; text-align: center; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
.stat-val { font-size: 1.6rem; font-weight: bold; color: #2563eb; }
.stat-label { color: #64748b; font-size: 0.75rem; text-transform: uppercase; letter-spacing: 0.05em; margin-top: 4px; }
table { width: 100%; border-collapse: collapse; margin: 1rem 0; font-size: 0.85rem; }
th, td { padding: 10px 12px; text-align: left; border-bottom: 1px solid #e2e8f0; }
th { background: #f1f5f9; font-weight: 600; color: #475569; cursor: pointer; }
tr:hover { background: #f8fafc; }
.plot-grid { display: grid; grid-template-columns: 1fr; gap: 2rem; margin: 2rem 0; }
.plot-grid img { max-width: 100%; height: auto; border-radius: 8px; border: 1px solid #e2e8f0; }
.plot-2col { display: grid; grid-template-columns: 1fr 1fr; gap: 2rem; margin: 2rem 0; }
.plot-2col img { max-width: 100%; height: auto; border-radius: 8px; border: 1px solid #e2e8f0; }
.section { margin: 2rem 0; }
.group-tag { display: inline-block; padding: 2px 8px; border-radius: 4px; font-size: 0.75rem; font-weight: 600; color: white; }
.footer { text-align: center; margin-top: 3rem; padding-top: 1.5rem; border-top: 1px solid #e2e8f0; color: #94a3b8; font-size: 0.8rem; }
.highlight-row { background-color: #fef3c7 !important; }
@media print {
    body { background: white; font-size: 10pt; }
    .container { box-shadow: none; padding: 0; }
    .header { -webkit-print-color-adjust: exact; }
    .plot-2col { display: block; }
    .plot-2col > div { margin-bottom: 2rem; }
}
"""

js = """
function sortTable(th) {
    const table = th.closest('table');
    const tbody = table.querySelector('tbody');
    const rows = Array.from(tbody.querySelectorAll('tr'));
    const idx = Array.from(th.parentNode.children).indexOf(th);
    const asc = th.dataset.sort !== 'asc';
    rows.sort((a, b) => {
        let va = a.children[idx].textContent.replace(/,/g, '');
        let vb = b.children[idx].textContent.replace(/,/g, '');
        const na = parseFloat(va), nb = parseFloat(vb);
        if (!isNaN(na) && !isNaN(nb)) return asc ? na - nb : nb - na;
        return asc ? va.localeCompare(vb) : vb.localeCompare(va);
    });
    rows.forEach(r => tbody.appendChild(r));
    th.dataset.sort = asc ? 'asc' : 'desc';
}
"""

# KPI numbers
n_samples = len(df_black) // 2
n_as = len(df[df['condition'] == 'AS'])
n_n = len(df[df['condition'] == 'N'])
mean_enrichment = df_enrichment['absolute_enrichment'].replace([np.inf], np.nan).mean() if len(df_enrichment) > 0 else 0
mean_rel_enrichment = df_enrichment['relative_enrichment'].replace([np.inf], np.nan).mean() if len(df_enrichment) > 0 else 0

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Black Sample Comparison Report</title>
    <style>{css}</style>
</head>
<body>
<div class="container">
    <div class="header">
        <h1>Adaptive Sampling Comparison Report</h1>
        <p>Black Samples — Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
        <div class="summary-stats" style="margin-top:1.5rem;">
            <div class="stat-card"><div class="stat-val">{n_samples}</div><div class="stat-label">Samples</div></div>
            <div class="stat-card"><div class="stat-val">{n_samples * 2}</div><div class="stat-label">Sequencing Runs</div></div>
            <div class="stat-card"><div class="stat-val">3</div><div class="stat-label">Extraction Methods</div></div>
            <div class="stat-card"><div class="stat-val">{mean_enrichment:.1f}×</div><div class="stat-label">Mean Enrichment</div></div>
        </div>
    </div>
"""

# --- Section 1: Experimental Design ---
html += """
<div class="section">
<h2>Experimental Design</h2>
<table>
<thead><tr><th>Group</th><th>Samples</th><th>Swab Type</th><th>Extraction Kit</th><th>Barcodes</th></tr></thead>
<tbody>
<tr><td><span class="group-tag" style="background:#3b82f6">A</span> Sponge / PowerSoil</td><td>A1–A5</td><td>Sponge</td><td>PowerSoil</td><td>3, 4, 5, 6, 7</td></tr>
<tr><td><span class="group-tag" style="background:#f59e0b">C</span> Cotton / Zymo</td><td>C1–C5</td><td>Cotton</td><td>Zymo</td><td>14, 15, 16, 18, 19</td></tr>
<tr><td><span class="group-tag" style="background:#8b5cf6">D</span> Zymo / Zymo</td><td>D1–D5</td><td>Zymo</td><td>Zymo</td><td>26, 27, 28, 29, 30</td></tr>
</tbody>
</table>
<p>Each sample was sequenced twice: once with <strong>Adaptive Sampling</strong> (AS, Nanopore real-time target enrichment) and once under <strong>Normal</strong> (N, standard) conditions.</p>
</div>
"""

# --- Section 2: QC Overview ---
html += '<div class="section"><h2>Sequencing Quality Overview</h2>'

# Build QC table
qc_display_cols = ['sample_id', 'group', 'condition', 'dna_concentration_ng_ul']
if reads_col:
    qc_display_cols.append(reads_col)
bases_col = None
for c in ['number_of_bases', 'Total Bases', 'total_bases']:
    if c in df.columns:
        bases_col = c
        break
if bases_col:
    qc_display_cols.append(bases_col)

existing_qc = [c for c in qc_display_cols if c in df.columns]
if existing_qc:
    qc_table = df[existing_qc].copy()
    qc_table = qc_table.rename(columns={
        'sample_id': 'Sample',
        'group': 'Group',
        'condition': 'Condition',
        'dna_concentration_ng_ul': 'DNA (ng/µL)',
    })
    if reads_col and reads_col in qc_table.columns:
        qc_table = qc_table.rename(columns={reads_col: 'Total Reads'})
    if bases_col and bases_col in qc_table.columns:
        qc_table = qc_table.rename(columns={bases_col: 'Total Bases'})

    qc_table = qc_table.sort_values(['Group', 'Sample', 'Condition'])
    html += df_to_html(qc_table, 'qc-table')

if 'total_reads' in plots:
    html += f'<div class="plot-grid"><img src="data:image/png;base64,{plots["total_reads"]}" alt="Total reads"></div>'
html += '</div>'

# --- Section 3: Adaptive Sampling Enrichment ---
html += '<div class="section"><h2>Adaptive Sampling Enrichment</h2>'
html += '<p>Enrichment ratio = Listeria reads (AS) / Listeria reads (N). Values &gt;1 indicate successful target enrichment.</p>'

if len(df_enrichment) > 0:
    enr_display = df_enrichment[['sample_id', 'group', 'swab_type', 'kit', 'dna_ng_ul',
                                  'listeria_AS', 'listeria_N', 'abs_difference',
                                  'pct_AS', 'pct_N',
                                  'absolute_enrichment', 'relative_enrichment',
                                  'as_mean_len', 'n_mean_len',
                                  'as_median_len', 'n_median_len']].copy()
    enr_display = enr_display.rename(columns={
        'sample_id': 'Sample', 'group': 'Group', 'swab_type': 'Swab',
        'kit': 'Kit', 'dna_ng_ul': 'DNA (ng/µL)',
        'listeria_AS': 'Listeria (AS)', 'listeria_N': 'Listeria (N)',
        'abs_difference': 'Δ Reads',
        'pct_AS': '% List. (AS)', 'pct_N': '% List. (N)',
        'absolute_enrichment': 'Abs. Enrich. (×)', 'relative_enrichment': 'Rel. Enrich. (×)',
        'as_mean_len': 'Mean Len (AS)', 'n_mean_len': 'Mean Len (N)',
        'as_median_len': 'Med. Len (AS)', 'n_median_len': 'Med. Len (N)',
    })
    html += df_to_html(enr_display, 'enrichment-table')

if 'abs_listeria' in plots:
    html += f'<div class="plot-grid"><img src="data:image/png;base64,{plots["abs_listeria"]}" alt="Absolute Listeria reads"></div>'
if 'pct_listeria' in plots:
    html += f'<div class="plot-grid"><img src="data:image/png;base64,{plots["pct_listeria"]}" alt="% Listeria"></div>'
if 'enrichment' in plots:
    html += f'<div class="plot-grid"><img src="data:image/png;base64,{plots["enrichment"]}" alt="Enrichment ratios"></div>'
if 'abs_vs_rel' in plots:
    html += f'<div class="plot-grid"><img src="data:image/png;base64,{plots["abs_vs_rel"]}" alt="Absolute vs Relative enrichment"></div>'
if 'group_enrichment' in plots:
    html += f'<div class="plot-grid"><img src="data:image/png;base64,{plots["group_enrichment"]}" alt="Group enrichment"></div>'
html += '</div>'

# --- Section 4: Extraction Method Comparison ---
html += '<div class="section"><h2>Extraction Method Comparison</h2>'
html += '<p>Comparison of sequencing yield and Listeria detection across swab types and DNA extraction kits.</p>'

# Summary stats per group
group_stats = []
for grp in ['Black_A', 'Black_C', 'Black_D']:
    sub = df[df['group'] == grp]
    meta_sub = df_black[df_black['group'] == grp].iloc[0]
    stat = {
        'Group': grp.replace('Black_', ''),
        'Swab': meta_sub['swab_type'],
        'Kit': meta_sub['kit'],
        'Mean DNA (ng/µL)': sub['dna_concentration_ng_ul'].mean(),
    }
    if reads_col and reads_col in sub.columns:
        stat['Mean Reads'] = int(sub[reads_col].mean())
    if listeria_col and listeria_col in sub.columns:
        stat['Mean Listeria'] = int(sub[listeria_col].mean())
    group_stats.append(stat)

df_group_stats = pd.DataFrame(group_stats)
html += df_to_html(df_group_stats, 'group-stats-table')

if 'extraction_comparison' in plots:
    html += f'<div class="plot-grid"><img src="data:image/png;base64,{plots["extraction_comparison"]}" alt="Extraction comparison"></div>'
if 'dna_vs_reads' in plots:
    html += f'<div class="plot-grid"><img src="data:image/png;base64,{plots["dna_vs_reads"]}" alt="DNA vs Reads"></div>'
html += '</div>'

# --- Section 5: AMR Overview ---
html += '<div class="section"><h2>AMR Gene Detection</h2>'

if len(df_amr_reads) > 0:
    # Filter to Black barcodes
    if 'Barcode' in df_amr_reads.columns:
        black_barcode_strs = [f'barcode{bc:02d}' for bc in BLACK_BARCODES]
        amr_black = df_amr_reads[df_amr_reads['Barcode'].isin(black_barcode_strs)].copy()
        if len(amr_black) > 0:
            html += '<h3>AMR Genes in Reads (Black Samples)</h3>'
            html += df_to_html(amr_black, 'amr-reads-table')
        else:
            html += '<p>No AMR genes detected in Black sample reads.</p>'
    else:
        html += '<p>AMR reads data format not recognized.</p>'
else:
    html += '<p>No AMR reads data available.</p>'

if len(df_amr_contigs) > 0:
    if 'Barcode' in df_amr_contigs.columns:
        amr_ctg_black = df_amr_contigs[df_amr_contigs['Barcode'].isin(black_barcode_strs)].copy()
        if len(amr_ctg_black) > 0:
            html += '<h3>AMR Genes in Contigs (Black Samples)</h3>'
            html += df_to_html(amr_ctg_black, 'amr-contigs-table')
html += '</div>'

# Footer
html += f"""
<div class="footer">
    <p>Comparison Report — Generated {datetime.now().strftime('%Y-%m-%d %H:%M')} — Nanopore Adaptive Sampling Pipeline</p>
</div>
</div>
<script>{js}</script>
</body>
</html>
"""

# Write output
report_path = os.path.join(out_dir, 'comparison_report.html')
with open(report_path, 'w') as f:
    f.write(html)
print(f"\nReport saved: {report_path}")
print("Done!")

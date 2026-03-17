#!/usr/bin/env python3
"""
Step 17 report builder (v2, overhauled).
Purpose: generate the main HTML report(s) used at the end of the pipeline.

Includes:
  - Round-aware report generation (per-round + combined)
  - Metadata-enriched tables (cohort, group, type)
  - Extraction method comparison section
  - Read statistics overview
  - Listeria analysis (reads + contigs)
  - Assembly statistics with Type column
  - AMR overview with Type/Cohort context
  - Interactive sortable tables
  - Embedded plots

Usage:
    python3 17_generate_report_v2.py <base_dir>
"""

import sys
import os
import re
import glob
import base64
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io

# Colors (matching local report palette)
C_AS = '#86efac'  # Green-300  (used in tables)
C_N  = '#fca5a5'  # Red-300    (used in tables)
C_NEUTRAL = '#e2e8f0'

# Condition colors (matching build_local_black_report.py)
COND_COLORS = {'AS': '#2d8a4e', 'N': '#b03060'}
COND_FILLS  = {'AS': '#2d8a4e55', 'N': '#b0306055'}

# Swab-family colors (matching local report group suffix colors)
SWAB_COLORS = {
    'Sponge':    '#3b82f6',  # blue  (suffix A)
    'Cotton':    '#f59e0b',  # amber (suffix C)
    'Zymo_swab': '#8b5cf6',  # purple (suffix D)
}

# Per swab/kit combo: shade varies by kit within each swab family
METHOD_COLORS = {
    'Sponge / PowerSoil': '#2563eb',
    'Sponge / Micro':     '#60a5fa',
    'Sponge / Mini':      '#1d4ed8',
    'Cotton / Zymo':      '#d97706',
    'Cotton / Micro':     '#fbbf24',
    'Cotton / Mini':      '#b45309',
    'Zymo_swab / Zymo':   '#7c3aed',
    'Zymo_swab / Micro':  '#a78bfa',
    'Zymo_swab / Mini':   '#5b21b6',
}
METHOD_FALLBACK_COLORS = ['#64748b', '#0f766e', '#b45309', '#be123c']

# ============================================================
# Plot helpers
# ============================================================

def create_boxplot(df, x_col, y_col, title, ylabel):
    """Generate base64 encoded boxplot comparing groups."""
    plt.figure(figsize=(8, 5))

    unique_vals = sorted(df[x_col].unique())
    if set(unique_vals) == {'AS', 'N'} or set(unique_vals) == {'N', 'AS'}:
        unique_vals = ['AS', 'N']

    data = [df[df[x_col] == v][y_col].dropna().tolist() for v in unique_vals]
    colors = [C_AS if 'AS' in str(v) else C_N for v in unique_vals]

    box = plt.boxplot(data, labels=unique_vals, patch_artist=True,
                      medianprops=dict(color='black', linewidth=1.5),
                      flierprops=dict(marker='o', markersize=4, alpha=0.5))

    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)

    plt.title(title)
    plt.ylabel(ylabel)
    plt.grid(False)
    plt.tight_layout()

    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight', dpi=100)
    plt.close()
    return base64.b64encode(buf.getvalue()).decode('utf-8')


def create_vertical_bar_chart(df, x_col, y_col, title, ylabel, ylim=None):
    """Generate bar chart of sums per group."""
    plt.figure(figsize=(8, 5))
    sums = df.groupby(x_col)[y_col].sum()
    labels = sorted(sums.index)
    if set(labels) == {'AS', 'N'} or set(labels) == {'N', 'AS'}:
        labels = ['AS', 'N']

    values = [sums[l] for l in labels]
    colors = [C_AS if 'AS' in str(l) else C_N for l in labels]

    plt.bar(labels, values, color=colors, alpha=0.9, edgecolor='black')
    plt.title(title)
    plt.ylabel(ylabel)
    if ylim:
        plt.ylim(0, ylim)
    plt.grid(False)

    for i, v in enumerate(values):
        plt.text(i, v, f'{int(v):,}', ha='center', va='bottom', fontsize=10)

    plt.tight_layout()
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight', dpi=100)
    plt.close()
    return base64.b64encode(buf.getvalue()).decode('utf-8')


def create_read_length_histogram(df_rl, title):
    """Generate a log-scale read length histogram, AS vs N overlaid. Returns base64 PNG."""
    if df_rl is None or len(df_rl) == 0:
        return None

    # Derive condition from sample name
    if 'condition' not in df_rl.columns:
        df_rl = df_rl.copy()
        df_rl['condition'] = df_rl['sample'].apply(extract_type)

    max_len = max(int(df_rl['length'].max()), 1000)
    bins = np.logspace(np.log10(50), np.log10(max(max_len, 100000)), 101)

    fig, ax = plt.subplots(figsize=(9, 4.5))
    for cond, color, edge in [('N', COND_FILLS['N'], COND_COLORS['N']), ('AS', COND_FILLS['AS'], COND_COLORS['AS'])]:
        sub = df_rl[df_rl['condition'] == cond]
        if len(sub) == 0:
            continue
        ax.hist(sub['length'], bins=bins, weights=sub['count'],
                histtype='stepfilled', color=color, edgecolor=edge,
                linewidth=0.9, alpha=0.55, label=cond)

    ax.set_xscale('log')
    ax.set_xlabel('Read length (bp)')
    ax.set_ylabel('Read count')
    ax.set_title(title)
    ax.set_xticks([100, 1000, 10000, 100000])
    ax.set_xticklabels(['100', '1k', '10k', '100k'])
    ax.axvline(400, linestyle='--', color='#6b7280', linewidth=1.0)
    ax.legend(frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=120)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode('utf-8')


def method_label_for_sample(meta_lookup, sample):
    """Return 'swab_type / kit' label for a single sample."""
    info = meta_lookup.get(sample, {})
    swab = info.get('swab_type', '') or ''
    kit = info.get('kit', '') or ''
    parts = [p for p in [swab, kit] if p and p not in ('None', 'NA', 'nan')]
    return ' / '.join(parts) if parts else ''


def get_method_color(method_label):
    """Get color for a swab/kit method combo."""
    if method_label in METHOD_COLORS:
        return METHOD_COLORS[method_label]
    # Fallback for unknown combos
    idx = hash(method_label) % len(METHOD_FALLBACK_COLORS)
    return METHOD_FALLBACK_COLORS[idx]


def swab_type_from_label(label):
    """Extract swab type from a 'Swab / Kit' label."""
    return label.split(' / ')[0] if ' / ' in label else label


def create_method_comparison_chart(method_df, group_col, type_col='type'):
    """Generate 3-panel grouped bar chart: total reads, listeria reads (log), listeria %.
    Returns base64 PNG or None."""
    if method_df is None or len(method_df) == 0:
        return None

    groups = sorted(method_df[group_col].dropna().unique())
    if len(groups) == 0:
        return None

    # Aggregate per group+type
    agg_cols = {}
    if 'number_of_reads' in method_df.columns:
        agg_cols['total_reads'] = ('number_of_reads', 'sum')
    if 'l_reads' in method_df.columns:
        agg_cols['listeria_reads'] = ('l_reads', 'sum')

    if not agg_cols:
        return None

    summary = method_df.groupby([group_col, type_col], as_index=False).agg(**agg_cols)

    if 'total_reads' in summary.columns and 'listeria_reads' in summary.columns:
        summary['listeria_pct'] = np.where(
            summary['total_reads'] > 0,
            summary['listeria_reads'] / summary['total_reads'] * 100.0,
            0.0
        )
    else:
        return None

    x = np.arange(len(groups))
    width = 0.36

    fig, axes = plt.subplots(1, 3, figsize=(15.0, 4.4), sharex=True)
    specs = [
        ('total_reads', 'Total filtered reads', False),
        ('listeria_reads', 'Listeria-classified reads', True),
        ('listeria_pct', 'Listeria fraction (%)', False),
    ]

    for ax, (metric, title, log_scale) in zip(axes, specs):
        n_vals, as_vals = [], []
        for group in groups:
            sub = summary[summary[group_col] == group]
            n_vals.append(float(sub.loc[sub[type_col] == 'N', metric].sum()))
            as_vals.append(float(sub.loc[sub[type_col] == 'AS', metric].sum()))

        ax.bar(x - width / 2, n_vals, width, color=COND_FILLS['N'], edgecolor=COND_COLORS['N'], label='N')
        ax.bar(x + width / 2, as_vals, width, color=COND_FILLS['AS'], edgecolor=COND_COLORS['AS'], label='AS')
        ax.set_title(title)
        ax.set_xticks(x)
        ax.set_xticklabels(groups, rotation=18, ha='right')
        if log_scale:
            ax.set_yscale('log')
        ax.spines['top'].set_visible(False)
        ax.spines['right'].set_visible(False)

    axes[0].legend(frameon=False)
    plt.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=120)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode('utf-8')


def create_method_boxplot(df, group_col, y_col, hue_col, title, ylabel):
    """Generate grouped boxplot: groups on x, AS/N side-by-side."""
    groups = sorted(df[group_col].dropna().unique())
    hues = ['AS', 'N']
    n_groups = len(groups)
    if n_groups == 0:
        return None

    fig, ax = plt.subplots(figsize=(max(8, n_groups * 2), 5))
    width = 0.35
    positions_as = np.arange(n_groups) - width / 2
    positions_n = np.arange(n_groups) + width / 2

    for hue, pos in [('AS', positions_as), ('N', positions_n)]:
        face = COND_FILLS[hue]
        edge = COND_COLORS[hue]
        data_list = []
        for g in groups:
            subset = df[(df[group_col] == g) & (df[hue_col] == hue)][y_col].dropna().tolist()
            data_list.append(subset if subset else [0])
        bp = ax.boxplot(data_list, positions=pos, widths=width * 0.8,
                        patch_artist=True, showfliers=False,
                        medianprops=dict(color=edge, linewidth=1.6),
                        boxprops=dict(facecolor=face, edgecolor=edge, linewidth=1.2),
                        whiskerprops=dict(color=edge, linewidth=1.0),
                        capprops=dict(color=edge, linewidth=1.0))

    ax.set_xticks(np.arange(n_groups))
    ax.set_xticklabels(groups, rotation=30, ha='right')
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend([plt.Rectangle((0, 0), 1, 1, fc=COND_FILLS['AS'], ec=COND_COLORS['AS']),
               plt.Rectangle((0, 0), 1, 1, fc=COND_FILLS['N'], ec=COND_COLORS['N'])],
              ['AS', 'N'], frameon=False)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()

    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=120)
    plt.close(fig)
    return base64.b64encode(buf.getvalue()).decode('utf-8')


# ============================================================
# Utility helpers
# ============================================================

def img_to_base64(path):
    """Embed image as base64 for self-contained HTML."""
    if not os.path.exists(path):
        return None
    with open(path, 'rb') as f:
        data = base64.b64encode(f.read()).decode('utf-8')
    ext = path.rsplit('.', 1)[-1]
    mime = 'image/png' if ext == 'png' else 'image/jpeg'
    return f'data:{mime};base64,{data}'


def df_to_html_table(df, table_id='', max_rows=None):
    """Convert DataFrame to nice HTML table."""
    if max_rows and len(df) > max_rows:
        df = df.head(max_rows)
    classes = 'data-table sortable'
    html = f'<table class="{classes}" id="{table_id}">\n<thead>\n<tr>'
    for col in df.columns:
        html += f'<th onclick="sortTable(this)">{col}</th>'
    html += '</tr>\n</thead>\n<tbody>\n'
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


def extract_type(sample_name):
    """Extract AS or N from sample name."""
    m = re.search(r'barcode\d+_(AS|N)', str(sample_name))
    return m.group(1) if m else ('AS' if '_AS' in str(sample_name) else 'N')


def extract_round(sample_name):
    """Extract round number from sample name (e.g., r1_barcode03_AS -> r1)."""
    m = re.match(r'^(r\d+)_', str(sample_name))
    return m.group(1) if m else None


# ============================================================
# Main
# ============================================================

base_dir = sys.argv[1]
out_dir = os.path.join(base_dir, 'processing/report')
os.makedirs(out_dir, exist_ok=True)

# ============================================================
# Load data
# ============================================================

print("Loading data...")

# Read metrics
read_csv = os.path.join(base_dir, 'processing/stats/read_metrics_summary.csv')
df_reads = pd.read_csv(read_csv) if os.path.exists(read_csv) else pd.DataFrame()

# Derive type column
if len(df_reads) > 0 and 'type' not in df_reads.columns and 'sample' in df_reads.columns:
    df_reads['type'] = df_reads['sample'].apply(extract_type)

# ---- Read length distributions ----
read_lengths_path = os.path.join(base_dir, 'processing/read_lengths_filtered_agg.tsv')
df_read_lengths = None
if os.path.exists(read_lengths_path):
    try:
        df_read_lengths = pd.read_csv(read_lengths_path, sep='\t', header=None,
                                       names=['sample', 'length', 'state', 'count'])
        # Keep only 'filtered' state if present
        if 'state' in df_read_lengths.columns:
            states = df_read_lengths['state'].unique()
            if 'filtered' in states:
                df_read_lengths = df_read_lengths[df_read_lengths['state'] == 'filtered'].copy()
        print(f"Loaded read length distributions: {len(df_read_lengths)} rows, "
              f"{df_read_lengths['sample'].nunique()} samples")
    except Exception as e:
        print(f"Warning: Failed to load read_lengths_filtered_agg.tsv: {e}")
        df_read_lengths = None
else:
    print(f"Note: No read length distribution file at {read_lengths_path}")

# ---- Metadata loading ----
df_meta = None
meta_path_1 = os.path.join(base_dir, 'sample_metadata.csv')
meta_path_2 = os.path.join(base_dir, 'samplesheets/sample_sheet_master.csv')

for mp in [meta_path_1, meta_path_2]:
    if os.path.exists(mp):
        try:
            df_meta = pd.read_csv(mp)
            print(f"Loaded metadata from: {mp}")
            # Ensure basename column exists for join
            if 'basename' not in df_meta.columns:
                # Try to derive it
                if all(c in df_meta.columns for c in ['round', 'barcode', 'condition']):
                    df_meta['basename'] = df_meta.apply(
                        lambda r: f"r{int(r['round'])}_barcode{int(r['barcode']):02d}_{r['condition']}",
                        axis=1
                    )
            break
        except Exception as e:
            print(f"Warning: Failed to load metadata from {mp}: {e}")
            df_meta = None

# Build a metadata lookup keyed by basename
meta_lookup = {}
if df_meta is not None and 'basename' in df_meta.columns:
    meta_cols = ['basename', 'cohort', 'group', 'swab_type', 'kit', 'round', 'condition']
    existing_meta_cols = [c for c in meta_cols if c in df_meta.columns]
    meta_sub = df_meta[existing_meta_cols].drop_duplicates(subset='basename')
    meta_lookup = meta_sub.set_index('basename').to_dict('index')
    print(f"Metadata lookup built with {len(meta_lookup)} entries.")


def enrich_with_metadata(df, sample_col='sample'):
    """Add Cohort, Group columns from metadata if available."""
    if df is None or len(df) == 0 or not meta_lookup:
        return df
    df = df.copy()
    if 'Cohort' not in df.columns:
        df['Cohort'] = df[sample_col].map(lambda s: meta_lookup.get(s, {}).get('cohort', ''))
    if 'Group' not in df.columns:
        df['Group'] = df[sample_col].map(lambda s: meta_lookup.get(s, {}).get('group', ''))
    return df


# ---- Listeria overview ----
listeria_csv = os.path.join(base_dir, 'processing/listeria/overview/listeria_overview.csv')
listeria_tsv = os.path.join(base_dir, 'processing/listeria/listeria_summary.tsv')

df_listeria = None

# 1. Try loading compiled CSV (Step 15 output)
if os.path.exists(listeria_csv):
    try:
        df_listeria = pd.read_csv(listeria_csv)
        print(f"Loaded listeria_overview.csv with columns: {df_listeria.columns.tolist()}")

        rename_map = {
            'Sample': 'sample',
            'Listeria Reads': 'listeria_reads',
            'Listeria (%)': 'listeria_ratio',
            'Mean Read Length': 'listeria_mean_len',
            'Median Read Length': 'listeria_median_len',
            'Type': 'type',
            'pct_listeria': 'listeria_ratio',
            'mean_read_len': 'listeria_mean_len',
            'median_read_len': 'listeria_median_len',
        }
        df_listeria.rename(columns=rename_map, inplace=True)

        if 'type' not in df_listeria.columns and 'sample' in df_listeria.columns:
            df_listeria['type'] = df_listeria['sample'].apply(extract_type)

    except Exception as e:
        print(f"Warning: Failed to load listeria_overview.csv: {e}")
        df_listeria = None

# 2. Fallback to raw TSV
if df_listeria is None and os.path.exists(listeria_tsv):
    try:
        print("Loading listeria_summary.tsv (fallback)...")
        df_listeria = pd.read_csv(listeria_tsv, sep='\t', header=None)

        if len(df_listeria.columns) >= 4 and isinstance(df_listeria.iloc[0, 1], (int, float, np.number)):
            cols = ['sample', 'listeria_reads', 'listeria_bases', 'listeria_mean_len']
            if len(df_listeria.columns) == 5:
                cols.append('listeria_median_len')
            df_listeria.columns = cols[:len(df_listeria.columns)]

        if 'sample' in df_listeria.columns:
            df_listeria['type'] = df_listeria['sample'].apply(extract_type)

        if 'listeria_ratio' not in df_listeria.columns:
            if len(df_reads) > 0 and 'sample' in df_reads.columns:
                temp = df_reads[['sample', 'number_of_reads']].copy()
                df_listeria = pd.merge(df_listeria, temp, on='sample', how='left')
                df_listeria['listeria_ratio'] = (
                    df_listeria['listeria_reads'] / df_listeria['number_of_reads'] * 100
                ).fillna(0)
            else:
                df_listeria['listeria_ratio'] = 0.0

    except Exception as e:
        print(f"Warning: Failed to load listeria summary TSV: {e}")
        df_listeria = None

# Final check and merge contigs
if df_listeria is not None:
    contigs_summary_tsv = os.path.join(base_dir, 'processing/listeria/listeria_contigs_summary.tsv')
    if os.path.exists(contigs_summary_tsv):
        try:
            df_lc = pd.read_csv(contigs_summary_tsv, sep='\t', header=None,
                                names=['sample', 'assembler', 'contigs', 'bases', 'median', 'k_count'])
            df_lc = df_lc.groupby(['sample', 'assembler'], as_index=False).max()

            df_lc_p = df_lc.pivot(index='sample', columns='assembler',
                                  values=['contigs', 'bases', 'median'])
            df_lc_p.columns = [f"{val}_{asm}" for val, asm in df_lc_p.columns]
            df_lc_p = df_lc_p.reset_index()
            df_lc_p.columns = [c.replace('metaMDBG', 'mdbg') for c in df_lc_p.columns]

            df_listeria = pd.merge(df_listeria, df_lc_p, on='sample', how='left')
            for c in df_lc_p.columns:
                if c != 'sample':
                    df_listeria[c] = df_listeria[c].fillna(0)
                    if 'contigs' in c or 'bases' in c:
                        df_listeria[c] = df_listeria[c].astype(int)
        except Exception as e:
            print(f"Warning: Failed to load listeria contigs summary: {e}")

# ---- AMR ----
amr_reads_csv = os.path.join(base_dir, 'processing/amrfinder/overview/amr_reads_overview.csv')
df_amr_reads = pd.read_csv(amr_reads_csv) if os.path.exists(amr_reads_csv) else None

amr_contigs_csv = os.path.join(base_dir, 'processing/amrfinder/overview/amr_contigs_overview.csv')
df_amr_contigs = pd.read_csv(amr_contigs_csv) if os.path.exists(amr_contigs_csv) else None

# Add type/cohort to AMR tables
for amr_df_name in ['df_amr_reads', 'df_amr_contigs']:
    amr_df = locals()[amr_df_name]
    if amr_df is not None:
        sample_col = 'Sample' if 'Sample' in amr_df.columns else ('sample' if 'sample' in amr_df.columns else None)
        if sample_col:
            if 'type' not in amr_df.columns and 'Type' not in amr_df.columns:
                amr_df['Type'] = amr_df[sample_col].apply(extract_type)
            if meta_lookup:
                if 'Cohort' not in amr_df.columns:
                    amr_df['Cohort'] = amr_df[sample_col].map(
                        lambda s: meta_lookup.get(s, {}).get('cohort', ''))
        locals()[amr_df_name] = amr_df

# ---- Assembly Stats ----
stats_mdbg_tsv = os.path.join(base_dir, 'processing/stats/assembly_stats_mdbg.tsv')
stats_flye_tsv = os.path.join(base_dir, 'processing/stats/assembly_stats_flye.tsv')
stats_myloasm_tsv = os.path.join(base_dir, 'processing/stats/assembly_stats_myloasm.tsv')

df_stats_mdbg = pd.read_csv(stats_mdbg_tsv, sep='\t') if os.path.exists(stats_mdbg_tsv) else None
df_stats_flye = pd.read_csv(stats_flye_tsv, sep='\t') if os.path.exists(stats_flye_tsv) else None
df_stats_myloasm = pd.read_csv(stats_myloasm_tsv, sep='\t') if os.path.exists(stats_myloasm_tsv) else None


def clean_assembly_stats(df_s):
    """Clean assembly stats: extract sample, add Type and Cohort columns."""
    if df_s is not None and 'file' in df_s.columns:
        df_s = df_s.copy()
        df_s['Sample'] = df_s['file'].apply(
            lambda x: os.path.basename(os.path.dirname(x)) if '/' in x else x
        )
        df_s['Sample'] = df_s['Sample'].str.replace('assembly_', '')
        df_s['Type'] = df_s['Sample'].apply(extract_type)

        if meta_lookup:
            df_s['Cohort'] = df_s['Sample'].map(
                lambda s: meta_lookup.get(s, {}).get('cohort', ''))

        cols_map = {'num_seqs': 'Contigs', 'N50': 'N50', 'Q2': 'Median Length'}
        base_keep = ['Sample', 'Type']
        if 'Cohort' in df_s.columns:
            base_keep.append('Cohort')
        keep = base_keep + list(cols_map.keys())
        keep = [c for c in keep if c in df_s.columns]
        return df_s[keep].rename(columns=cols_map)
    return df_s


if df_stats_mdbg is not None:
    df_stats_mdbg = clean_assembly_stats(df_stats_mdbg)
if df_stats_flye is not None:
    df_stats_flye = clean_assembly_stats(df_stats_flye)
if df_stats_myloasm is not None:
    df_stats_myloasm = clean_assembly_stats(df_stats_myloasm)

# ---- Listeria plots (pre-generated PNGs) ----
overview_dir = os.path.join(base_dir, 'processing/listeria/overview')
plots = {}
for name in ['pct_listeria_per_barcode', 'listeria_reads_log_per_barcode',
             'listeria_contigs_comparison', 'listeria_multi_panel']:
    png = os.path.join(overview_dir, f'{name}.png')
    b64 = img_to_base64(png)
    if b64:
        plots[name] = b64

# ============================================================
# Round detection
# ============================================================

all_samples = set()
if len(df_reads) > 0 and 'sample' in df_reads.columns:
    all_samples.update(df_reads['sample'].tolist())
if df_listeria is not None and 'sample' in df_listeria.columns:
    all_samples.update(df_listeria['sample'].tolist())

detected_rounds = set()
for s in all_samples:
    r = extract_round(s)
    if r:
        detected_rounds.add(r)

detected_rounds = sorted(detected_rounds)
print(f"Detected rounds: {detected_rounds if detected_rounds else 'none (single round)'}")


def filter_by_round(df, round_prefix, sample_col='sample'):
    """Filter dataframe to samples matching a given round prefix."""
    if df is None or len(df) == 0:
        return df
    if sample_col not in df.columns:
        return df
    mask = df[sample_col].str.startswith(f'{round_prefix}_')
    return df[mask].copy()


# ============================================================
# PDF Generation
# ============================================================

def generate_pdf(html_path):
    """Generate PDF using headless Chrome if available."""
    import subprocess
    pdf_path = html_path.replace('.html', '.pdf')
    chrome_paths = [
        '/Applications/Google Chrome.app/Contents/MacOS/Google Chrome',
        '/usr/bin/google-chrome',
        '/usr/bin/chromium-browser'
    ]

    chrome_bin = None
    for path in chrome_paths:
        if os.path.exists(path):
            chrome_bin = path
            break

    if chrome_bin:
        print(f"Generating PDF: {pdf_path} (using Chrome)...")
        try:
            cmd = f'"{chrome_bin}" --headless --disable-gpu --print-to-pdf="{pdf_path}" "{html_path}"'
            subprocess.run(cmd, shell=True, check=True, stderr=subprocess.DEVNULL)
            print(f"PDF created successfully: {pdf_path}")
        except Exception as e:
            print(f"Warning: PDF generation failed: {e}")
    else:
        print("Note: PDF generation skipped (Chrome/Chromium not found).")


# ============================================================
# CSS
# ============================================================

css_style = """
body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 0; background-color: #f8fafc; color: #334155; }
.container { max-width: 1200px; margin: 0 auto; padding: 20px; background-color: #ffffff; box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1); }
h1, h2, h3 { color: #1e293b; border-bottom: 2px solid #e2e8f0; padding-bottom: 10px; margin-top: 2rem; }
.header { text-align: center; margin-bottom: 2rem; padding: 2rem; background: linear-gradient(135deg, #3b82f6 0%, #2563eb 100%); color: white; border-radius: 8px; }
.header h1 { color: white; border: none; margin: 0; }
.summary-stats { display: flex; gap: 20px; flex-wrap: wrap; margin-bottom: 2rem; justify-content: center; }
.stat-card { flex: 1; min-width: 150px; max-width: 200px; background: #fff; border: 1px solid #e2e8f0; border-radius: 8px; padding: 1rem; text-align: center; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
.stat-val { font-size: 1.8rem; font-weight: bold; color: #2563eb; line-height: 1.2; }
.stat-label { color: #64748b; font-size: 0.8rem; text-transform: uppercase; letter-spacing: 0.05em; margin-top: 5px; }
table { width: 100%; border-collapse: collapse; margin: 1rem 0; font-size: 0.9rem; }
th, td { padding: 12px 15px; text-align: left; border-bottom: 1px solid #e2e8f0; }
th { background-color: #f1f5f9; font-weight: 600; color: #475569; cursor: pointer; }
tr:hover { background-color: #f8fafc; }
.plot-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(450px, 1fr)); gap: 2rem; margin: 2rem 0; }
.plot-grid img { max-width: 100%; height: auto; border-radius: 8px; border: 1px solid #e2e8f0; }
.footer { text-align: center; margin-top: 4rem; padding-top: 2rem; border-top: 1px solid #e2e8f0; color: #94a3b8; font-size: 0.85rem; }
.tabs { display: flex; border-bottom: 1px solid #e2e8f0; margin-bottom: 1rem; }
.tab-btn { padding: 10px 20px; background: none; border: none; border-bottom: 2px solid transparent; cursor: pointer; font-weight: 600; color: #64748b; }
.tab-btn:hover { color: #2563eb; }
.tab-btn.active { border-bottom-color: #2563eb; color: #2563eb; }
.tab-content { display: none; animation: fadeEffect 0.5s; }
.tab-content.active { display: block; }
@keyframes fadeEffect { from {opacity: 0;} to {opacity: 1;} }
.search-box { margin-bottom: 1rem; }
.search-box input { padding: 8px 12px; border: 1px solid #e2e8f0; border-radius: 4px; width: 300px; }
.table-wrapper { overflow-x: auto; max-height: 500px; }
p.note { color: #64748b; font-size: 0.9rem; margin: 0.5rem 0 1.5rem 0; line-height: 1.5; }

/* Table of contents - pill links */
.toc { display: flex; flex-wrap: wrap; gap: 10px; justify-content: center; margin: 1.5rem 0 2rem 0; padding: 1rem; }
.toc a { display: inline-block; padding: 8px 18px; background: #eff6ff; color: #2563eb; border-radius: 999px; text-decoration: none; font-size: 0.85rem; font-weight: 600; border: 1px solid #bfdbfe; transition: background 0.2s; }
.toc a:hover { background: #dbeafe; }

/* Print / PDF Styles */
@media print {
    body { background: white; font-size: 10pt; }
    .container { max-width: 100%; box-shadow: none; padding: 0; }
    .header { background: #2563eb !important; -webkit-print-color-adjust: exact; color: white !important; }
    .header h1 { color: white !important; }
    .stat-card { border: 1px solid #ccc; box-shadow: none; break-inside: avoid; }
    .tab-content { display: block !important; opacity: 1 !important; margin-bottom: 2rem; }
    .tabs, .tab-btn { display: none !important; }
    .search-box { display: none !important; }
    .toc { display: none !important; }
    .footer { margin-top: 2rem; border-top: 1px solid #ccc; }
    .table-wrapper { max-height: none !important; overflow: visible !important; }
    table { page-break-inside: auto; }
    tr { page-break-inside: avoid; page-break-after: auto; }
    thead { display: table-header-group; }
    h2 { page-break-before: auto; margin-top: 2rem; border-bottom: 1px solid #000; }
    .section { margin-bottom: 2rem; break-inside: avoid; }
    .plot-grid { display: block; }
    .plot-grid > div { margin-bottom: 2rem; break-inside: avoid; }
    .plot-grid img { max-width: 80%; margin: 0 auto; display: block; }
}
"""

# ============================================================
# JavaScript
# ============================================================

js_code = """
function openTab(evt, tabName) {
  var i, tabcontent, tablinks;
  // Find the parent section of the clicked tab
  var parentSection = evt.currentTarget.closest('.section');
  if (!parentSection) parentSection = document;
  tabcontent = parentSection.querySelectorAll(".tab-content");
  for (i = 0; i < tabcontent.length; i++) {
    tabcontent[i].style.display = "none";
    tabcontent[i].classList.remove("active");
  }
  tablinks = parentSection.querySelectorAll(".tab-btn");
  for (i = 0; i < tablinks.length; i++) {
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }
  document.getElementById(tabName).style.display = "block";
  document.getElementById(tabName).classList.add("active");
  evt.currentTarget.className += " active";
}

function sortTable(th) {
  var table = th.closest("table");
  var tbody = table.querySelector("tbody");
  var rows = Array.from(tbody.querySelectorAll("tr"));
  var index = Array.from(th.parentNode.children).indexOf(th);
  var asc = th.getAttribute("data-order") !== "asc";

  rows.sort(function(a, b) {
    var valA = a.children[index].innerText;
    var valB = b.children[index].innerText;
    var numA = parseFloat(valA.replace(/,/g, ""));
    var numB = parseFloat(valB.replace(/,/g, ""));
    if (!isNaN(numA) && !isNaN(numB)) {
      return asc ? numA - numB : numB - numA;
    }
    return asc ? valA.localeCompare(valB) : valB.localeCompare(valA);
  });

  rows.forEach(function(row) { tbody.appendChild(row); });
  th.setAttribute("data-order", asc ? "asc" : "desc");
}

function filterTable(input, tableId) {
  var filter = input.value.toUpperCase();
  var table = document.getElementById(tableId);
  var tr = table.getElementsByTagName("tr");
  for (var i = 1; i < tr.length; i++) {
    var td = tr[i].getElementsByTagName("td")[0];
    if (td) {
      var txtValue = td.textContent || td.innerText;
      if (txtValue.toUpperCase().indexOf(filter) > -1) {
        tr[i].style.display = "";
      } else {
        tr[i].style.display = "none";
      }
    }
  }
}
"""


# ============================================================
# Report builder function
# ============================================================

def build_report(df_reads_r, df_listeria_r, df_amr_reads_r, df_amr_contigs_r,
                 df_stats_mdbg_r, df_stats_flye_r, df_stats_myloasm_r,
                 round_label=None, df_read_lengths_r=None):
    """Build a single HTML report string for the given data subset.

    round_label: e.g. "Round 1" or None for combined.
    """

    # --- Summary stats ---
    n_total = len(df_reads_r) if df_reads_r is not None and len(df_reads_r) > 0 else 0
    n_as = 0
    n_n = 0
    total_reads = 0
    total_bases_gb = 0.0

    if df_reads_r is not None and len(df_reads_r) > 0:
        if 'type' in df_reads_r.columns:
            n_as = len(df_reads_r[df_reads_r['type'] == 'AS'])
            n_n = len(df_reads_r[df_reads_r['type'] == 'N'])
        if 'number_of_reads' in df_reads_r.columns:
            total_reads = df_reads_r['number_of_reads'].sum()
        if 'total_bases' in df_reads_r.columns:
            total_bases_gb = df_reads_r['total_bases'].sum() / 1e9

    listeria_samples = 0
    pos_as = 0
    pos_n = 0
    if df_listeria_r is not None and 'listeria_reads' in df_listeria_r.columns:
        positives = df_listeria_r[df_listeria_r['listeria_reads'] > 0]
        listeria_samples = len(positives)
        if 'type' in positives.columns:
            pos_as = len(positives[positives['type'] == 'AS'])
            pos_n = len(positives[positives['type'] == 'N'])

    genes_set = set()
    if df_amr_reads_r is not None and 'Gene Symbol' in df_amr_reads_r.columns:
        genes_set.update(df_amr_reads_r['Gene Symbol'].unique())
    if df_amr_contigs_r is not None and 'Gene Symbol' in df_amr_contigs_r.columns:
        genes_set.update(df_amr_contigs_r['Gene Symbol'].unique())
    n_genes = len(genes_set)

    s_total_reads = f"{total_reads/1e6:.1f}M" if total_reads > 1e6 else f"{total_reads:,}"
    s_total_bases = f"{total_bases_gb:.2f}Gb"
    s_pos_breakdown = f"{listeria_samples} ({pos_as} AS, {pos_n} N)"

    title_suffix = f" - {round_label}" if round_label else ""
    report_title = f"Listeria Analysis Report{title_suffix}"

    # --- HTML head & header ---
    html_out = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>{report_title}</title>
    <style>{css_style}</style>
</head>
<body>
<div class="container">
    <div class="header">
        <h1>{report_title}</h1>
        <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>

        <div class="summary-stats" style="margin-top: 2rem;">
            <div class="stat-card">
                <div class="stat-val">{n_total}</div>
                <div class="stat-label">Total Samples</div>
            </div>
            <div class="stat-card">
                <div class="stat-val">{s_total_reads}</div>
                <div class="stat-label">Total Reads</div>
            </div>
            <div class="stat-card">
                <div class="stat-val">{s_total_bases}</div>
                <div class="stat-label">Total Bases</div>
            </div>
            <div class="stat-card">
                <div class="stat-val">{s_pos_breakdown}</div>
                <div class="stat-label">Listeria Positive</div>
            </div>
            <div class="stat-card">
                <div class="stat-val">{n_genes}</div>
                <div class="stat-label">Resistance Genes</div>
            </div>
        </div>
    </div>
"""

    # --- Table of Contents ---
    html_out += """
    <div class="toc">
        <a href="#sec-qc">QC Overview</a>
        <a href="#sec-read-lengths">Read Lengths</a>
        <a href="#sec-listeria-reads">Listeria Reads</a>
        <a href="#sec-method-comparison">Method Comparison</a>
        <a href="#sec-listeria-contigs">Listeria Contigs</a>
        <a href="#sec-assembly">Assembly Statistics</a>
        <a href="#sec-amr">AMR Detection</a>
    </div>
"""

    html_body = ""

    # ================================================================
    # 1. QC Section
    # ================================================================
    html_body += '<div class="section" id="sec-qc"><h2>Sequencing Quality Control</h2>'
    html_body += '<p class="note">Read-level statistics after adapter trimming and length filtering. Each row represents one demultiplexed barcode library.</p>'
    html_body += '<p class="note">Use this section first to judge scale before looking at per-barcode or per-method differences.</p>'
    if df_reads_r is not None and len(df_reads_r) > 0:
        disp_qc = df_reads_r.copy()
        # Enrich with metadata
        if 'sample' in disp_qc.columns:
            disp_qc = enrich_with_metadata(disp_qc, 'sample')

        cols_qc = ['sample', 'type']
        if 'Cohort' in disp_qc.columns:
            cols_qc.append('Cohort')
        cols_qc += ['number_of_reads', 'mean_read_length', 'median_read_length',
                     'read_length_N50', 'total_bases']
        existing_qc = [c for c in cols_qc if c in disp_qc.columns]
        disp_qc = disp_qc[existing_qc].rename(columns={
            'sample': 'Sample', 'type': 'Type',
            'number_of_reads': 'Total Reads',
            'mean_read_length': 'Mean Length (bp)',
            'median_read_length': 'Median Length (bp)',
            'read_length_N50': 'Read N50 (bp)',
            'total_bases': 'Total Bases (bp)'
        })

        for col in ['Mean Length (bp)', 'Median Length (bp)']:
            if col in disp_qc.columns:
                disp_qc[col] = disp_qc[col].round(1)

        html_body += '<div class="search-box"><input type="text" placeholder="Search samples..." '
        html_body += 'onkeyup="filterTable(this, \'qc-table\')"></div>'
        html_body += f'<div class="table-wrapper">{df_to_html_table(disp_qc, "qc-table")}</div>'
    else:
        html_body += "<p>No read statistics available.</p>"
    html_body += '</div>'

    # ================================================================
    # 1b. Read Length Distributions
    # ================================================================
    html_body += '<div class="section" id="sec-read-lengths"><h2>Read Length Distributions</h2>'
    html_body += '<p class="note">All barcodes are pooled here so the AS and N read-length profiles can be compared directly. The dashed line marks 400 bp (minimum length filter). Log-scale x-axis.</p>'
    html_body += '<p class="note">If one condition has a visibly larger area under the curve, it contributed more reads in that length range.</p>'

    if df_read_lengths_r is not None and len(df_read_lengths_r) > 0:
        # Overall AS vs N histogram
        rl_b64 = create_read_length_histogram(
            df_read_lengths_r,
            'Filtered read length distribution (all samples pooled)')
        if rl_b64:
            html_body += '<div style="max-width:900px; margin:1rem auto 2rem auto;">'
            html_body += f'<img src="data:image/png;base64,{rl_b64}" style="width:100%; border:1px solid #e2e8f0; border-radius:8px;">'
            html_body += '</div>'

        # Per-method read length comparison (if metadata available)
        if meta_lookup:
            rl_with_group = df_read_lengths_r.copy()
            if 'condition' not in rl_with_group.columns:
                rl_with_group['condition'] = rl_with_group['sample'].apply(extract_type)
            rl_with_group['method_label'] = rl_with_group['sample'].map(
                lambda s: method_label_for_sample(meta_lookup, s))
            rl_with_group = rl_with_group[rl_with_group['method_label'] != ''].copy()

            if len(rl_with_group) > 0:
                methods = sorted(rl_with_group['method_label'].unique())
                max_len = max(int(rl_with_group['length'].max()), 1000)
                bins = np.logspace(np.log10(50), np.log10(max(max_len, 100000)), 101)

                # --- Plot 1: 2-panel overview (N | AS), all methods overlaid ---
                # Matches section_plot_method_comparison from local report
                fig, axes = plt.subplots(1, 2, figsize=(12.8, 4.4), sharex=True, sharey=True)
                for ax, cond in zip(axes, ['N', 'AS']):
                    for method in methods:
                        sub = rl_with_group[
                            (rl_with_group['method_label'] == method) &
                            (rl_with_group['condition'] == cond)]
                        if len(sub) > 0:
                            color = get_method_color(method)
                            ax.hist(sub['length'], bins=bins, weights=sub['count'],
                                    histtype='stepfilled', color=color, edgecolor=color,
                                    linewidth=1.0, alpha=0.28, label=method)
                    ax.set_xscale('log')
                    ax.set_xlabel('Read length (bp)')
                    ax.set_title('Normal' if cond == 'N' else 'Adaptive Sampling')
                    ax.set_xticks([100, 1000, 10000, 100000])
                    ax.set_xticklabels(['100', '1k', '10k', '100k'])
                    ax.axvline(400, linestyle='--', color='#6b7280', linewidth=1.0)
                    ax.spines['top'].set_visible(False)
                    ax.spines['right'].set_visible(False)
                axes[0].set_ylabel('Read count')
                axes[0].legend(frameon=False)
                fig.suptitle('Read Length Comparison by Extraction Method', fontsize=12, y=1.02)
                buf = io.BytesIO()
                fig.savefig(buf, format='png', bbox_inches='tight', dpi=120)
                plt.close(fig)
                overview_b64 = base64.b64encode(buf.getvalue()).decode('utf-8')

                html_body += '<h3>Read Length Comparison by Extraction Method</h3>'
                html_body += '<p class="note">Two pooled comparison panels split by sequencing mode. Within each panel, all swab/kit methods are overlaid directly. Read left-to-right: first compare methods within Normal, then within Adaptive Sampling.</p>'
                html_body += '<div style="max-width:100%; margin:1rem auto 2rem auto;">'
                html_body += f'<img src="data:image/png;base64,{overview_b64}" style="width:100%; border:1px solid #e2e8f0; border-radius:8px;">'
                html_body += '</div>'

                # --- Plot 2: per swab-type panels, kits overlaid, AS/N with condition fills ---
                # Matches section_plot_read_lengths_by_method style from local report
                swab_types = sorted(set(swab_type_from_label(m) for m in methods))
                n_panels = len(swab_types)
                if n_panels > 0:
                    fig, axes = plt.subplots(1, n_panels, figsize=(5.0 * max(n_panels, 1), 4.2),
                                             sharey=True)
                    axes = np.atleast_1d(axes).ravel()

                    for ax, swab in zip(axes, swab_types):
                        swab_methods = sorted(m for m in methods
                                              if swab_type_from_label(m) == swab)
                        for method in swab_methods:
                            kit = method.split(' / ')[1] if ' / ' in method else method
                            sub_method = rl_with_group[rl_with_group['method_label'] == method]
                            for cond in ['N', 'AS']:
                                sub = sub_method[sub_method['condition'] == cond]
                                if len(sub) > 0:
                                    ax.hist(sub['length'], bins=bins, weights=sub['count'],
                                            histtype='stepfilled',
                                            color=COND_FILLS[cond], edgecolor=COND_COLORS[cond],
                                            linewidth=1.0, alpha=0.5,
                                            label=f'{kit} ({cond})')
                        ax.set_xscale('log')
                        ax.set_xlabel('Read length (bp)')
                        ax.set_title(swab)
                        ax.set_xticks([100, 1000, 10000, 100000])
                        ax.set_xticklabels(['100', '1k', '10k', '100k'])
                        ax.axvline(400, linestyle='--', color='#6b7280', linewidth=1.0)
                        ax.spines['top'].set_visible(False)
                        ax.spines['right'].set_visible(False)
                    axes[0].set_ylabel('Read count')
                    axes[0].legend(frameon=False, fontsize=7)
                    fig.suptitle('Read Lengths by Swab Type (Kit Breakdown)', fontsize=12, y=1.02)
                    buf = io.BytesIO()
                    fig.savefig(buf, format='png', bbox_inches='tight', dpi=120)
                    plt.close(fig)
                    detail_b64 = base64.b64encode(buf.getvalue()).decode('utf-8')

                    html_body += '<h3>Read Lengths by Swab Type</h3>'
                    html_body += '<p class="note">Per swab type, kits overlaid with N and AS on the same axes. This view shows whether one extraction workflow tends to produce a different read-length profile overall.</p>'
                    html_body += '<div style="max-width:100%; margin:1rem auto 2rem auto;">'
                    html_body += f'<img src="data:image/png;base64,{detail_b64}" style="width:100%; border:1px solid #e2e8f0; border-radius:8px;">'
                    html_body += '</div>'
    else:
        html_body += '<p>No read length distribution data available (step 3b).</p>'

    html_body += '</div>'

    # ================================================================
    # 2. Listeria Reads Analysis
    # ================================================================
    html_body += '<div class="section" id="sec-listeria-reads"><h2>Listeria Reads Analysis</h2>'
    html_body += '<p class="note"><em>Listeria monocytogenes</em> reads identified by Kraken2 taxonomic classification. Percentage is relative to total filtered reads per sample.</p>'
    html_body += '<p class="note">This section is about enrichment visibility, not absolute sequencing yield. A small highlighted fraction means <em>Listeria</em> makes up only a small part of that sample even if the total library is large.</p>'
    html_body += '<p class="note">If you need the exact counts and percentages, use the sortable table directly below the figure.</p>'

    # Plots
    html_body += '<h3>Listeria Reads Comparison (AS vs N)</h3><div class="plot-grid">'
    comp_plots = []

    if df_listeria_r is not None and len(df_listeria_r) > 0:
        # Positive samples count
        pos_counts_df = pd.DataFrame([
            {'type': 'AS', 'count': pos_as},
            {'type': 'N', 'count': pos_n}
        ])
        p0 = create_vertical_bar_chart(pos_counts_df, 'type', 'count',
                                       'Listeria Positive Samples (AS vs N)', 'Count')
        comp_plots.append(('Positive Samples Count', p0))

    if df_reads_r is not None and len(df_reads_r) > 0 and 'type' in df_reads_r.columns:
        if 'median_read_length' in df_reads_r.columns:
            p1 = create_boxplot(df_reads_r, 'type', 'median_read_length',
                                'Median Read Length (AS vs N)', 'Length (bp)')
            comp_plots.append(('Median Read Length', p1))

    if df_listeria_r is not None and len(df_listeria_r) > 0:
        # Dynamic ylim instead of hardcoded 15000
        max_val = df_listeria_r.groupby('type')['listeria_reads'].sum().max()
        dynamic_ylim = int(max_val * 1.2) if max_val > 0 else None
        p3 = create_vertical_bar_chart(df_listeria_r, 'type', 'listeria_reads',
                                       'Total Listeria Reads (AS vs N)', 'Count',
                                       ylim=dynamic_ylim)
        comp_plots.append(('Total Listeria Reads', p3))

    for title, b64 in comp_plots:
        html_body += f'<div><h4>{title}</h4><img src="data:image/png;base64,{b64}"></div>'
    html_body += '</div>'

    # Listeria Reads Table
    if df_listeria_r is not None and len(df_listeria_r) > 0:
        html_body += '<h3>Per-Sample Listeria Reads Data</h3>'
        disp_listeria = df_listeria_r.copy()
        if 'sample' in disp_listeria.columns:
            disp_listeria = enrich_with_metadata(disp_listeria, 'sample')

        display_cols = ['sample', 'type']
        if 'Cohort' in disp_listeria.columns:
            display_cols.append('Cohort')
        display_cols += ['listeria_ratio', 'listeria_reads', 'listeria_bases',
                         'listeria_mean_len', 'listeria_median_len']
        existing_cols = [c for c in display_cols if c in disp_listeria.columns]

        final_df = disp_listeria[existing_cols].rename(columns={
            'sample': 'Sample', 'type': 'Type',
            'listeria_reads': 'Listeria Reads',
            'listeria_bases': 'Listeria Bases',
            'listeria_ratio': '% of Total Reads',
            'listeria_mean_len': 'Mean Length',
            'listeria_median_len': 'Median Length'
        })

        if '% of Total Reads' in final_df.columns:
            final_df['% of Total Reads'] = final_df['% of Total Reads'].map('{:.2f}'.format)

        html_body += '<div class="search-box"><input type="text" placeholder="Search samples..." '
        html_body += 'onkeyup="filterTable(this, \'listeria-reads-table\')"></div>'
        html_body += f'<div class="table-wrapper">{df_to_html_table(final_df, "listeria-reads-table")}</div>'
    html_body += '</div>'

    # ================================================================
    # 3. Extraction Method Comparison (NEW)
    # ================================================================
    html_body += '<div class="section" id="sec-method-comparison"><h2>Extraction Method Comparison</h2>'
    html_body += '<p class="note">Direct comparison of extraction methods (swab type + kit combinations). This section collapses barcode-level noise and shows how the methods behave as groups.</p>'
    html_body += '<p class="note">Use it to ask whether one extraction workflow systematically gives more total reads, more <em>Listeria</em> reads, or a higher enrichment fraction for Adaptive Sampling.</p>'

    has_method_data = (meta_lookup and df_listeria_r is not None and len(df_listeria_r) > 0
                       and df_reads_r is not None and len(df_reads_r) > 0)

    if has_method_data:
        # Build method-level summary
        method_df = df_reads_r.copy()
        if 'sample' in method_df.columns:
            method_df['Swab'] = method_df['sample'].map(
                lambda s: meta_lookup.get(s, {}).get('swab_type', ''))
            method_df['Kit'] = method_df['sample'].map(
                lambda s: meta_lookup.get(s, {}).get('kit', ''))
            method_df['method_label'] = method_df['sample'].map(
                lambda s: method_label_for_sample(meta_lookup, s))

        # Merge listeria reads
        if 'sample' in method_df.columns and df_listeria_r is not None:
            listeria_merge = df_listeria_r[['sample', 'listeria_reads']].copy()
            listeria_merge = listeria_merge.rename(columns={'listeria_reads': 'l_reads'})
            method_df = pd.merge(method_df, listeria_merge, on='sample', how='left')
            method_df['l_reads'] = method_df['l_reads'].fillna(0).astype(int)

        # Filter to rows that have method info
        method_df = method_df[method_df['method_label'] != ''].copy()

        if len(method_df) > 0:
            # Method-level totals table (group by swab+kit, not cohort-specific Group)
            method_summary = method_df.groupby(['Swab', 'Kit', 'type']).agg(
                Libraries=('sample', 'count'),
                Total_Reads=('number_of_reads', 'sum') if 'number_of_reads' in method_df.columns else ('sample', 'count'),
                Listeria_Reads=('l_reads', 'sum')
            ).reset_index()

            if 'Total_Reads' in method_summary.columns and 'Listeria_Reads' in method_summary.columns:
                method_summary['% Listeria'] = (
                    method_summary['Listeria_Reads'] / method_summary['Total_Reads'] * 100
                ).fillna(0).round(3)

            method_summary = method_summary.rename(columns={
                'type': 'Type', 'Total_Reads': 'Total Reads',
                'Listeria_Reads': 'Listeria Reads'
            })

            html_body += '<h3>Method-Level Summary</h3>'
            html_body += f'<div class="table-wrapper">{df_to_html_table(method_summary, "method-table")}</div>'

            # 3-panel grouped bar chart: total reads, listeria reads (log), listeria %
            if 'l_reads' in method_df.columns and 'method_label' in method_df.columns:
                chart_b64 = create_method_comparison_chart(method_df, 'method_label', 'type')
                if chart_b64:
                    html_body += '<h3>Extraction Method Comparison (3-Panel)</h3>'
                    html_body += '<div style="max-width:100%; margin:1rem auto 2rem auto;">'
                    html_body += f'<img src="data:image/png;base64,{chart_b64}" style="width:100%; border:1px solid #e2e8f0; border-radius:8px;">'
                    html_body += '</div>'

                # Boxplot: AS vs N listeria reads per method
                bp_b64 = create_method_boxplot(
                    method_df, 'method_label', 'l_reads', 'type',
                    'Listeria Read Counts by Method (AS vs N)',
                    'Listeria Reads'
                )
                if bp_b64:
                    html_body += '<h3>Listeria Reads by Extraction Method</h3>'
                    html_body += f'<div style="max-width:800px; margin:1rem auto 2rem auto;">'
                    html_body += f'<img src="data:image/png;base64,{bp_b64}" style="width:100%; border:1px solid #e2e8f0; border-radius:8px;">'
                    html_body += '</div>'
        else:
            html_body += '<p>No method-level metadata available for the current samples.</p>'
    else:
        html_body += '<p>No metadata file found. Method comparison requires sample_metadata.csv.</p>'
    html_body += '</div>'

    # ================================================================
    # 4. Listeria Contigs
    # ================================================================
    html_body += '<div class="section" id="sec-listeria-contigs"><h2>Listeria Contigs</h2>'
    html_body += '<p class="note">Contigs classified as <em>Listeria monocytogenes</em> from assembled genomes. Assemblers: Flye (metagenomic mode), MetaMDBG, Myloasm.</p>'
    html_body += '<p class="note">This section shifts from reads to assemblies. Use it to check whether the assemblies are centered on <em>Listeria</em> or still dominated by background organisms.</p>'

    has_contig_cols = (df_listeria_r is not None and len(df_listeria_r) > 0 and
                       any(c in df_listeria_r.columns for c in
                           ['contigs_flye', 'contigs_mdbg', 'contigs_myloasm']))

    if has_contig_cols:
        display_cols_contigs = ['sample', 'type',
                                'contigs_flye', 'bases_flye', 'median_flye',
                                'contigs_mdbg', 'bases_mdbg', 'median_mdbg',
                                'contigs_myloasm', 'bases_myloasm', 'median_myloasm']
        existing_cols_contigs = [c for c in display_cols_contigs if c in df_listeria_r.columns]

        final_df_contigs = df_listeria_r[existing_cols_contigs].rename(columns={
            'sample': 'Sample', 'type': 'Type',
            'contigs_flye': 'Flye Contigs', 'bases_flye': 'Flye Bases',
            'median_flye': 'Flye Median',
            'contigs_mdbg': 'MDBG Contigs', 'bases_mdbg': 'MDBG Bases',
            'median_mdbg': 'MDBG Median',
            'contigs_myloasm': 'Myloasm Contigs', 'bases_myloasm': 'Myloasm Bases',
            'median_myloasm': 'Myloasm Median'
        })

        # Enrich with Cohort
        if meta_lookup and 'Sample' in final_df_contigs.columns:
            final_df_contigs.insert(
                2, 'Cohort',
                final_df_contigs['Sample'].map(
                    lambda s: meta_lookup.get(s, {}).get('cohort', ''))
            )

        # Filter out rows where all contig counts are 0
        contig_count_cols = [c for c in ['Flye Contigs', 'MDBG Contigs', 'Myloasm Contigs']
                             if c in final_df_contigs.columns]
        if contig_count_cols:
            final_df_contigs = final_df_contigs[
                final_df_contigs[contig_count_cols].sum(axis=1) > 0
            ].copy()

        html_body += '<div class="search-box"><input type="text" placeholder="Search samples..." '
        html_body += 'onkeyup="filterTable(this, \'listeria-contigs-table\')"></div>'
        html_body += f'<div class="table-wrapper">{df_to_html_table(final_df_contigs, "listeria-contigs-table")}</div>'
    else:
        html_body += "<p>No Listeria contig counts available (run Step 14).</p>"
    html_body += '</div>'

    # ================================================================
    # 5. Assembly Statistics
    # ================================================================
    html_body += '<div class="section" id="sec-assembly"><h2>Assembly Statistics</h2>'
    html_body += '<p class="note">Genome assembly metrics for all contigs (not filtered to <em>Listeria</em>). The assembly-metric panels use a log10 y-axis so differences between lower and higher values stay visible.</p>'
    html_body += '<p class="note">Zero-valued assemblies are omitted from log-scale boxplots. Type column derived from sample name (AS = adaptive sampling, N = normal).</p>'
    html_body += """
  <div class="tabs">
    <button class="tab-btn active" onclick="openTab(event, 'stats-mdbg')">MetaMDBG</button>
    <button class="tab-btn" onclick="openTab(event, 'stats-flye')">Flye</button>
    <button class="tab-btn" onclick="openTab(event, 'stats-myloasm')">Myloasm</button>
  </div>
"""

    for tab_id, tab_label, df_s in [
        ('stats-mdbg', 'MetaMDBG', df_stats_mdbg_r),
        ('stats-flye', 'Flye', df_stats_flye_r),
        ('stats-myloasm', 'Myloasm', df_stats_myloasm_r),
    ]:
        active = ' active' if tab_id == 'stats-mdbg' else ''
        html_body += f'<div class="tab-content{active}" id="{tab_id}">'
        if df_s is not None and len(df_s) > 0:
            html_body += f'<div class="table-wrapper">{df_to_html_table(df_s, f"{tab_id}-table")}</div>'
        else:
            html_body += f"<p>No {tab_label} stats available.</p>"
        html_body += '</div>'

    html_body += '</div>'

    # ================================================================
    # 6. AMR Gene Detection
    # ================================================================
    html_body += '<div class="section" id="sec-amr"><h2>AMR Gene Detection</h2>'
    html_body += '<p class="note">Antimicrobial resistance genes identified by AMRFinderPlus. Detected in both raw reads (FASTA) and assembled contigs.</p>'
    html_body += '<p class="note">The AMR tables keep barcode-level rows so you can see which sample each call came from, rather than only showing cohort-wide totals.</p>'
    html_body += """
  <div class="tabs">
    <button class="tab-btn active" onclick="openTab(event, 'amr-reads')">AMR Genes in Reads</button>
    <button class="tab-btn" onclick="openTab(event, 'amr-contigs')">AMR Genes in Contigs</button>
  </div>
"""

    # AMR from Reads
    html_body += '<div class="tab-content active" id="amr-reads">'
    if df_amr_reads_r is not None and len(df_amr_reads_r) > 0:
        sample_col_amr = 'Sample' if 'Sample' in df_amr_reads_r.columns else 'sample'
        type_col_amr = 'Type' if 'Type' in df_amr_reads_r.columns else 'type'

        if 'Gene Symbol' in df_amr_reads_r.columns and type_col_amr in df_amr_reads_r.columns:
            amr_counts = df_amr_reads_r.groupby(type_col_amr)['Gene Symbol'].nunique().reset_index()
            amr_counts.columns = ['type', 'Gene Count']
            p_amr = create_vertical_bar_chart(amr_counts, 'type', 'Gene Count',
                                              'Unique AMR Genes in Reads (AS vs N)', 'Count')
            html_body += f'<div style="max-width:600px; margin-bottom:2rem;">'
            html_body += f'<img src="data:image/png;base64,{p_amr}" style="width:100%; border:1px solid #e2e8f0; border-radius:4px;">'
            html_body += '</div>'

        html_body += f'<div class="table-wrapper">{df_to_html_table(df_amr_reads_r, "amr-reads-table", max_rows=500)}</div>'
    else:
        html_body += "<p>No AMR genes found in reads.</p>"
    html_body += '</div>'

    # AMR from Contigs
    html_body += '<div class="tab-content" id="amr-contigs">'
    if df_amr_contigs_r is not None and len(df_amr_contigs_r) > 0:
        sample_col_amr = 'Sample' if 'Sample' in df_amr_contigs_r.columns else 'sample'
        type_col_amr = 'Type' if 'Type' in df_amr_contigs_r.columns else 'type'

        if 'Gene Symbol' in df_amr_contigs_r.columns and type_col_amr in df_amr_contigs_r.columns:
            amr_counts_c = df_amr_contigs_r.groupby(type_col_amr)['Gene Symbol'].nunique().reset_index()
            amr_counts_c.columns = ['type', 'Gene Count']
            p_amr_c = create_vertical_bar_chart(amr_counts_c, 'type', 'Gene Count',
                                                'Unique AMR Genes in Contigs (AS vs N)', 'Count')
            html_body += f'<div style="max-width:600px; margin-bottom:2rem;">'
            html_body += f'<img src="data:image/png;base64,{p_amr_c}" style="width:100%; border:1px solid #e2e8f0; border-radius:4px;">'
            html_body += '</div>'

        html_body += f'<div class="table-wrapper">{df_to_html_table(df_amr_contigs_r, "amr-contigs-table")}</div>'
    else:
        html_body += "<p>No AMR genes found in contigs.</p>"
    html_body += '</div></div>'

    # ================================================================
    # Footer
    # ================================================================
    html_footer = f"""
<div class="footer">
  Raw data processed by Tim Thilo Maria Reska's Pipeline v1.0.<br>
  Report generated via customized Python engine.
</div>

</div> <!-- End Container -->

<script>
{js_code}
</script>
</body>
</html>
"""

    return html_out + html_body + html_footer


# ============================================================
# Filter assembly stats by round
# ============================================================

def filter_assembly_stats_by_round(df_s, round_prefix):
    """Filter assembly stats df (already cleaned) to a specific round."""
    if df_s is None or len(df_s) == 0 or 'Sample' not in df_s.columns:
        return df_s
    mask = df_s['Sample'].str.startswith(f'{round_prefix}_')
    result = df_s[mask].copy()
    return result if len(result) > 0 else None


def filter_amr_by_round(df_amr, round_prefix):
    """Filter AMR df to samples in a given round."""
    if df_amr is None or len(df_amr) == 0:
        return df_amr
    sample_col = 'Sample' if 'Sample' in df_amr.columns else ('sample' if 'sample' in df_amr.columns else None)
    if sample_col is None:
        return df_amr
    mask = df_amr[sample_col].str.startswith(f'{round_prefix}_')
    result = df_amr[mask].copy()
    return result if len(result) > 0 else None


def filter_by_cohort(df, cohort, sample_col='sample'):
    """Filter dataframe to samples belonging to a specific cohort via meta_lookup."""
    if df is None or len(df) == 0 or not meta_lookup:
        return df
    if sample_col not in df.columns:
        return df
    mask = df[sample_col].map(lambda s: meta_lookup.get(s, {}).get('cohort', '')) == cohort
    result = df[mask].copy()
    return result if len(result) > 0 else None


def filter_asm_by_cohort(df_s, cohort):
    """Filter cleaned assembly stats by cohort (uses Sample column)."""
    if df_s is None or len(df_s) == 0 or 'Sample' not in df_s.columns or not meta_lookup:
        return df_s
    mask = df_s['Sample'].map(lambda s: meta_lookup.get(s, {}).get('cohort', '')) == cohort
    result = df_s[mask].copy()
    return result if len(result) > 0 else None


def filter_amr_by_cohort(df_amr, cohort):
    """Filter AMR df to samples in a specific cohort."""
    if df_amr is None or len(df_amr) == 0 or not meta_lookup:
        return df_amr
    sample_col = 'Sample' if 'Sample' in df_amr.columns else ('sample' if 'sample' in df_amr.columns else None)
    if sample_col is None:
        return df_amr
    mask = df_amr[sample_col].map(lambda s: meta_lookup.get(s, {}).get('cohort', '')) == cohort
    result = df_amr[mask].copy()
    return result if len(result) > 0 else None


# ============================================================
# Detect cohorts
# ============================================================

detected_cohorts = set()
if meta_lookup:
    for info in meta_lookup.values():
        c = info.get('cohort', '')
        if c and c != 'Control':
            detected_cohorts.add(c)
detected_cohorts = sorted(detected_cohorts)
print(f"Detected cohorts: {detected_cohorts if detected_cohorts else 'none'}")


# ============================================================
# Generate reports
# ============================================================

print("Generating HTML report(s)...")

# --- 1. Combined report (all rounds, all cohorts) ---
combined_html = build_report(
    df_reads, df_listeria, df_amr_reads, df_amr_contigs,
    df_stats_mdbg, df_stats_flye, df_stats_myloasm,
    round_label=None, df_read_lengths_r=df_read_lengths
)

out_file = os.path.join(out_dir, 'pipeline_report.html')
with open(out_file, 'w') as f:
    f.write(combined_html)
print(f"Report generated: {out_file}")

# Export CSVs from combined data
if len(df_reads) > 0:
    qc_csv = os.path.join(out_dir, 'qc_metrics.csv')
    df_reads.to_csv(qc_csv, index=False)
    print(f"Exported QC table to: {qc_csv}")

if df_listeria is not None and len(df_listeria) > 0:
    listeria_reads_csv = os.path.join(out_dir, 'listeria_reads_summary.csv')
    df_listeria.to_csv(listeria_reads_csv, index=False)
    print(f"Exported Listeria Reads table to: {listeria_reads_csv}")

for label, df_s in [('mdbg', df_stats_mdbg), ('flye', df_stats_flye), ('myloasm', df_stats_myloasm)]:
    if df_s is not None and len(df_s) > 0:
        csv_path = os.path.join(out_dir, f'assembly_stats_{label}.csv')
        df_s.to_csv(csv_path, index=False)
        print(f"Exported {label} stats to: {csv_path}")

if df_amr_reads is not None:
    df_amr_reads.to_csv(os.path.join(out_dir, 'amr_genes_reads.csv'), index=False)
if df_amr_contigs is not None:
    df_amr_contigs.to_csv(os.path.join(out_dir, 'amr_genes_contigs.csv'), index=False)

# --- 2. Per-round reports (all cohorts combined, one per round) ---
if len(detected_rounds) > 1:
    for rnd in detected_rounds:
        round_num = rnd.replace('r', '')
        round_label = f"Round {round_num}"

        dr = filter_by_round(df_reads, rnd) if len(df_reads) > 0 else pd.DataFrame()
        dl = filter_by_round(df_listeria, rnd) if df_listeria is not None else None
        ar = filter_amr_by_round(df_amr_reads, rnd)
        ac = filter_amr_by_round(df_amr_contigs, rnd)
        sm = filter_assembly_stats_by_round(df_stats_mdbg, rnd)
        sf = filter_assembly_stats_by_round(df_stats_flye, rnd)
        sy = filter_assembly_stats_by_round(df_stats_myloasm, rnd)
        rl = filter_by_round(df_read_lengths, rnd) if df_read_lengths is not None else None

        round_html = build_report(dr, dl, ar, ac, sm, sf, sy, round_label=round_label,
                                  df_read_lengths_r=rl)
        round_file = os.path.join(out_dir, f'pipeline_report_{rnd}.html')
        with open(round_file, 'w') as f:
            f.write(round_html)
        print(f"Round report generated: {round_file}")

# --- 3. Per-cohort per-round reports (Black r1, Black r2, Blue r1, Blue r2) ---
if detected_cohorts:
    # Determine which round labels to iterate
    round_list = detected_rounds if len(detected_rounds) > 1 else [None]

    for cohort in detected_cohorts:
        for rnd in round_list:
            # Start from full data, filter by cohort
            dr_c = filter_by_cohort(df_reads, cohort) if len(df_reads) > 0 else pd.DataFrame()
            dl_c = filter_by_cohort(df_listeria, cohort) if df_listeria is not None else None
            ar_c = filter_amr_by_cohort(df_amr_reads, cohort)
            ac_c = filter_amr_by_cohort(df_amr_contigs, cohort)
            sm_c = filter_asm_by_cohort(df_stats_mdbg, cohort)
            sf_c = filter_asm_by_cohort(df_stats_flye, cohort)
            sy_c = filter_asm_by_cohort(df_stats_myloasm, cohort)
            rl_c = filter_by_cohort(df_read_lengths, cohort) if df_read_lengths is not None else None

            # Then filter by round if applicable
            if rnd is not None:
                round_num = rnd.replace('r', '')
                label = f"{cohort} - Round {round_num}"
                suffix = f"_{cohort.lower()}_{rnd}"

                dr_c = filter_by_round(dr_c, rnd) if dr_c is not None and len(dr_c) > 0 else pd.DataFrame()
                dl_c = filter_by_round(dl_c, rnd) if dl_c is not None else None
                ar_c = filter_amr_by_round(ar_c, rnd)
                ac_c = filter_amr_by_round(ac_c, rnd)
                sm_c = filter_assembly_stats_by_round(sm_c, rnd)
                sf_c = filter_assembly_stats_by_round(sf_c, rnd)
                sy_c = filter_assembly_stats_by_round(sy_c, rnd)
                rl_c = filter_by_round(rl_c, rnd) if rl_c is not None else None
            else:
                label = cohort
                suffix = f"_{cohort.lower()}"

            cohort_html = build_report(dr_c, dl_c, ar_c, ac_c, sm_c, sf_c, sy_c, round_label=label,
                                       df_read_lengths_r=rl_c)
            cohort_file = os.path.join(out_dir, f'pipeline_report{suffix}.html')
            with open(cohort_file, 'w') as f:
                f.write(cohort_html)
            print(f"Cohort report generated: {cohort_file}")

# --- 4. Generate PDFs ---
generate_pdf(out_file)
if len(detected_rounds) > 1:
    for rnd in detected_rounds:
        generate_pdf(os.path.join(out_dir, f'pipeline_report_{rnd}.html'))

print("Done.")

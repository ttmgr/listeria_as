#!/usr/bin/env python3
"""
Step 21 statistical analysis.
Purpose: run formal tests that compare extraction/swab methods and AS vs N outcomes.
Starts from sample metadata and extends automatically when comparison tables exist.

Usage:
    python3 21_statistical_analysis.py <base_dir>

Output:
    stats_results.txt  — plain text summary
    stats_report.html  — HTML report with tables, plots, and method scorecard
"""

import sys
import os
import pandas as pd
import numpy as np
from scipy import stats
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io
import base64
from itertools import combinations
from datetime import datetime

# ============================================================
# Config
# ============================================================

BLACK_BARCODES = [6, 7, 8, 9, 10, 16, 17, 18, 19, 20, 26, 27, 28, 29, 30, 34]
GROUP_LABELS = {
    'Black_A': 'A',
    'Black_C': 'C',
    'Black_D': 'D',
}
GROUP_COLORS = {'Black_A': '#3b82f6', 'Black_C': '#f59e0b', 'Black_D': '#8b5cf6'}

plt.rcParams.update({
    'font.family': 'sans-serif',
    'font.size': 11,
    'axes.titlesize': 13,
    'figure.dpi': 150,
})

# ============================================================
# Helpers
# ============================================================

def fig_to_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=150)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode('utf-8')


def kruskal_wallis(groups_dict, variable):
    """Run Kruskal-Wallis + pairwise Mann-Whitney on a dict of group→values."""
    results = {}
    arrays = [v for v in groups_dict.values() if len(v) > 0]
    labels = [k for k, v in groups_dict.items() if len(v) > 0]

    if len(arrays) < 2:
        return results

    # Overall KW test
    h, p_kw = stats.kruskal(*arrays)
    results['kruskal'] = {
        'H': round(h, 3),
        'p': round(p_kw, 4),
        'significant': p_kw < 0.05,
        'n_groups': len(arrays),
    }

    # Pairwise Mann-Whitney with Bonferroni correction
    pairs = list(combinations(range(len(arrays)), 2))
    n_comparisons = len(pairs)
    pairwise = []
    for i, j in pairs:
        u, p_mw = stats.mannwhitneyu(arrays[i], arrays[j], alternative='two-sided')
        p_corrected = min(p_mw * n_comparisons, 1.0)
        # Effect size: rank-biserial correlation
        n1, n2 = len(arrays[i]), len(arrays[j])
        r = 1 - (2 * u) / (n1 * n2)
        pairwise.append({
            'comparison': f'{labels[i]} vs {labels[j]}',
            'U': round(u, 1),
            'p_raw': round(p_mw, 4),
            'p_bonferroni': round(p_corrected, 4),
            'effect_r': round(r, 3),
            'significant': p_corrected < 0.05,
        })
    results['pairwise'] = pairwise
    return results


def wilcoxon_paired(as_vals, n_vals, label):
    """Wilcoxon signed-rank test for paired AS vs N."""
    if len(as_vals) < 4:
        return None
    diffs = np.array(as_vals) - np.array(n_vals)
    if np.all(diffs == 0):
        return {'label': label, 'W': 0, 'p': 1.0, 'significant': False,
                'median_AS': np.median(as_vals), 'median_N': np.median(n_vals),
                'note': 'All differences zero'}
    w, p = stats.wilcoxon(as_vals, n_vals)
    return {
        'label': label,
        'W': round(w, 1),
        'p': round(p, 4),
        'significant': p < 0.05,
        'median_AS': round(np.median(as_vals), 3),
        'median_N': round(np.median(n_vals), 3),
    }


def effect_size_d(a, b):
    """Cohen's d (pooled SD)."""
    na, nb = len(a), len(b)
    if na < 2 or nb < 2:
        return float('nan')
    pooled_sd = np.sqrt(((na - 1) * np.std(a, ddof=1)**2 + (nb - 1) * np.std(b, ddof=1)**2) / (na + nb - 2))
    if pooled_sd == 0:
        return float('nan')
    return round((np.mean(a) - np.mean(b)) / pooled_sd, 3)


# ============================================================
# Load data
# ============================================================

base_dir = sys.argv[1] if len(sys.argv) > 1 else '.'

# Try cluster path first, then preprocessing subdir
meta_path = None
for candidate in [
    os.path.join(base_dir, 'sample_metadata.csv'),
    os.path.join(base_dir, 'preprocessing', 'sample_metadata.csv'),
    os.path.join(base_dir, 'BAM', 'sample_metadata.csv'),
]:
    if os.path.exists(candidate):
        meta_path = candidate
        break

if not meta_path:
    print("ERROR: sample_metadata.csv not found. Run from project root or preprocessing dir.")
    sys.exit(1)

print(f"Loading metadata from: {meta_path}")
df_meta = pd.read_csv(meta_path)

# Filter to Black samples only (unique barcodes, deduplicate per barcode)
df_black = df_meta[df_meta['barcode'].isin(BLACK_BARCODES)].copy()
df_black['dna_concentration_ng_ul'] = pd.to_numeric(df_black['dna_concentration_ng_ul'], errors='coerce').fillna(0)
print(f"Black samples: {len(df_black)} rows ({df_black['barcode'].nunique()} barcodes × 2 conditions)")

for group_name in list(GROUP_LABELS):
    sub = df_black[df_black['group'] == group_name]
    if len(sub) == 0:
        continue
    parts = []
    if 'swab_type' in sub.columns:
        swabs = [str(v) for v in sub['swab_type'].dropna().unique() if str(v).strip()]
        if swabs:
            parts.append('/'.join(swabs))
    if 'kit' in sub.columns:
        kits = [str(v) for v in sub['kit'].dropna().unique() if str(v).strip()]
        if kits:
            parts.append('/'.join(kits))
    if parts:
        GROUP_LABELS[group_name] = f"{group_name.replace('Black_', '')} ({' / '.join(parts)})"

# One row per barcode (DNA concentration is the same for AS and N — same physical sample)
df_barcode = df_black.drop_duplicates(subset='barcode').copy()

# ============================================================
# Analysis 1: DNA concentration across groups
# ============================================================

print("\n--- Analysis 1: DNA concentration across groups ---")

groups_dna = {
    grp: df_barcode[df_barcode['group'] == grp]['dna_concentration_ng_ul'].values
    for grp in ['Black_A', 'Black_C', 'Black_D']
}
group_counts = [len(vals) for vals in groups_dna.values() if len(vals) > 0]
if group_counts and len(set(group_counts)) == 1:
    group_count_label = f"n={group_counts[0]} per group"
elif group_counts:
    group_count_label = "group sizes: " + ', '.join(str(v) for v in group_counts)
else:
    group_count_label = "group sizes unavailable"

# Descriptive stats
desc_rows = []
for grp, vals in groups_dna.items():
    # Exclude "lower than blank" zeros
    vals_nonzero = vals[vals > 0]
    desc_rows.append({
        'Group': GROUP_LABELS.get(grp, grp),
        'n': len(vals),
        'n (>0)': len(vals_nonzero),
        'Mean (ng/µL)': round(np.mean(vals), 3),
        'Median (ng/µL)': round(np.median(vals), 3),
        'SD': round(np.std(vals, ddof=1), 3) if len(vals) > 1 else 0,
        'Min': round(np.min(vals), 3),
        'Max': round(np.max(vals), 3),
    })
df_desc = pd.DataFrame(desc_rows)
print(df_desc.to_string(index=False))

# Kruskal-Wallis + pairwise
kw_results = kruskal_wallis(groups_dna, 'DNA concentration')
print(f"\nKruskal-Wallis: H={kw_results['kruskal']['H']}, p={kw_results['kruskal']['p']}", 
      "← SIGNIFICANT" if kw_results['kruskal']['significant'] else "← not significant")
print("\nPairwise comparisons (Bonferroni-corrected):")
df_pairwise = pd.DataFrame(kw_results['pairwise'])
print(df_pairwise.to_string(index=False))

# ============================================================
# Analysis 2: Spearman correlation DNA conc vs volume
# ============================================================

print("\n--- Analysis 2: Spearman DNA conc vs volume ---")
r, p = stats.spearmanr(df_barcode['dna_concentration_ng_ul'], df_barcode['volume_ul'])
print(f"Spearman r={round(r, 3)}, p={round(p, 4)}")

# Kit comparison: adapt to whatever the metadata actually contains.
print("\n--- Analysis 3: DNA concentration by extraction kit ---")
kit_groups = {
    str(kit): df_barcode[df_barcode['kit'] == kit]['dna_concentration_ng_ul'].values
    for kit in sorted(df_barcode['kit'].dropna().unique())
}
kit_rows = []
for kit_name, vals in kit_groups.items():
    kit_rows.append({
        'Kit': kit_name,
        'n': len(vals),
        'Median (ng/µL)': round(np.median(vals), 3) if len(vals) else np.nan,
        'Mean (ng/µL)': round(np.mean(vals), 3) if len(vals) else np.nan,
    })
df_kit_desc = pd.DataFrame(kit_rows)
print(df_kit_desc.to_string(index=False))

kit_test_label = 'Only one extraction kit present in metadata; no between-kit statistical test run.'
kit_significant = False
kit_test_type = 'not_run'
d = float('nan')

if len(kit_groups) == 2:
    kit_names = list(kit_groups)
    vals_a = kit_groups[kit_names[0]]
    vals_b = kit_groups[kit_names[1]]
    u, p_kit = stats.mannwhitneyu(vals_a, vals_b, alternative='two-sided')
    d = effect_size_d(vals_a, vals_b)
    kit_test_label = (
        f"Mann-Whitney U={round(u,1)}, p={round(p_kit, 4)}, Cohen's d={d}"
    )
    kit_significant = p_kit < 0.05
    kit_test_type = 'mannwhitney'
elif len(kit_groups) > 2:
    kw_kit = kruskal_wallis(kit_groups, 'DNA by kit')
    if 'kruskal' in kw_kit:
        kit_test_label = (
            f"Kruskal-Wallis H={kw_kit['kruskal']['H']}, p={kw_kit['kruskal']['p']}"
        )
        kit_significant = kw_kit['kruskal']['significant']
        kit_test_type = 'kruskal'
print(kit_test_label)

# Swab type comparison
print("\n--- Analysis 4: DNA concentration by swab type ---")
swab_groups = {
    swab: df_barcode[df_barcode['swab_type'] == swab]['dna_concentration_ng_ul'].values
    for swab in df_barcode['swab_type'].unique()
    if swab != '-'
}
kw_swab = kruskal_wallis(swab_groups, 'DNA by swab')
if 'kruskal' in kw_swab:
    print(f"Kruskal-Wallis: H={kw_swab['kruskal']['H']}, p={kw_swab['kruskal']['p']}",
          "← SIGNIFICANT" if kw_swab['kruskal']['significant'] else "← not significant")

# ============================================================
# METHOD SCORECARD — loads pipeline data when available
# ============================================================

# Metrics: (column_in_csv, display_name, higher_is_better)
SCORECARD_METRICS = [
    ('number_of_reads',        'Total Reads',            True),
    ('listeria_reads',         'Listeria Reads',          True),
    ('listeria_ratio',         '% Listeria / Total',      True),
    ('absolute_enrichment',    'Enrichment (Absolute)',   True),
    ('relative_enrichment',    'Enrichment (Relative)',   True),
    ('as_mean_len',            'Mean Read Length (AS)',   True),
    ('n_mean_len',             'Mean Read Length (N)',    True),
    ('as_median_len',          'Median Read Length (AS)', True),
    ('n_median_len',           'Median Read Length (N)',  True),
    ('Flye Contigs',           'Flye Contig Count',       True),
    ('MetaMDBG Contigs',       'MetaMDBG Contig Count',   True),
    ('Myloasm Contigs',        'Myloasm Contig Count',    True),
    ('flye_contig_bases',      'Flye Contig Bases',       True),
    ('mdbg_contig_bases',      'MetaMDBG Contig Bases',   True),
    ('myloasm_contig_bases',   'Myloasm Contig Bases',    True),
    ('flye_median_contig_len', 'Flye Median Contig Len',  True),
    ('mdbg_median_contig_len', 'MDBG Median Contig Len',  True),
    ('myloasm_median_contig_len','Myloasm Median Contig Len', True),
]

scorecard_html = ''
scorecard_rows = []
df_pipeline = None
df_enrich = None

# Try to find comparison_data.csv from pipeline
for cand in [
    os.path.join(base_dir, 'processing', 'report', 'comparison_data.csv'),
    os.path.join(base_dir, 'comparison_data.csv'),
]:
    if os.path.exists(cand):
        df_pipeline = pd.read_csv(cand)
        print(f"\nLoaded pipeline data: {cand} ({len(df_pipeline)} rows)")
        break

for cand in [
    os.path.join(base_dir, 'processing', 'report', 'enrichment_ratios.csv'),
    os.path.join(base_dir, 'enrichment_ratios.csv'),
]:
    if os.path.exists(cand):
        df_enrich = pd.read_csv(cand)
        print(f"Loaded enrichment data: {cand}")
        break

if df_pipeline is not None:
    print("\n--- Method Scorecard ---")
    # Merge enrichment data if available
    if df_enrich is not None and 'barcode' in df_enrich.columns and 'barcode' in df_pipeline.columns:
        merge_cols = [c for c in df_enrich.columns if c not in df_pipeline.columns or c == 'barcode']
        df_pipeline = df_pipeline.merge(df_enrich[merge_cols], on='barcode', how='left')

    # Filter Black samples only, per-condition
    if 'barcode' in df_pipeline.columns:
        df_pipeline = df_pipeline[df_pipeline['barcode'].isin(BLACK_BARCODES)].copy()

    # Standardize group column
    for gcol in ['group', 'Group']:
        if gcol in df_pipeline.columns:
            df_pipeline['_group'] = df_pipeline[gcol]
            break

    GROUPS = ['Black_A', 'Black_C', 'Black_D']
    GROUP_SHORT = {g: GROUP_LABELS.get(g, g.replace('Black_', '')) for g in GROUPS}
    scores = {g: 0 for g in GROUPS}

    for col, label, higher_better in SCORECARD_METRICS:
        if col not in df_pipeline.columns:
            continue
        df_pipeline[col] = pd.to_numeric(df_pipeline[col], errors='coerce')
        grp_vals = {
            g: df_pipeline[df_pipeline['_group'] == g][col].dropna().values
            for g in GROUPS
        }
        grp_medians = {g: np.median(v) if len(v) > 0 else np.nan for g, v in grp_vals.items()}
        valid = {g: m for g, m in grp_medians.items() if not np.isnan(m)}
        if not valid:
            continue

        # KW test
        arrays = [grp_vals[g] for g in GROUPS if len(grp_vals[g]) > 0]
        h_stat, p_val = (0, 1)
        if len(arrays) >= 2 and all(len(a) > 0 for a in arrays):
            try:
                h_stat, p_val = stats.kruskal(*arrays)
            except Exception:
                pass

        winner_grp = max(valid, key=lambda g: valid[g]) if higher_better else min(valid, key=lambda g: valid[g])
        scores[winner_grp] += 1

        row = {'Metric': label}
        for g in GROUPS:
            m = grp_medians.get(g, np.nan)
            row[GROUP_SHORT[g]] = f"{m:,.1f}" if not np.isnan(m) else '—'
        row['Winner'] = GROUP_SHORT[winner_grp]
        row['KW p'] = f"{p_val:.4f}"
        row['Significant'] = '✅' if p_val < 0.05 else '—'
        scorecard_rows.append(row)
        print(f"  {label:<35} winner={GROUP_SHORT[winner_grp]:<20} KW p={p_val:.4f}")

    # Total score row
    print(f"\nScores: {scores}")
    score_winner = max(scores, key=scores.get)
    print(f"Overall winner: {GROUP_SHORT[score_winner]} ({scores[score_winner]} metrics)")

    df_scorecard = pd.DataFrame(scorecard_rows)

    # Scorecard HTML
    scorecard_html = '<h2>📊 Method Scorecard</h2>'
    scorecard_html += ('<p>Median values per group for each metric. '
                       'Winner = group with highest (or best) value. '
                       'Highlighted rows are statistically significant (KW p&lt;0.05).</p>')

    # Score summary cards
    scorecard_html += '<div style="display:flex;gap:16px;margin:1rem 0;">'
    for g in GROUPS:
        color = GROUP_COLORS[g]
        scorecard_html += (f'<div style="flex:1;background:{color}22;border:2px solid {color};'
                           f'border-radius:8px;padding:1rem;text-align:center;">'
                           f'<div style="font-size:1.8rem;font-weight:bold;color:{color}">{scores[g]}</div>'
                           f'<div style="font-size:0.8rem;color:#475569;margin-top:4px">{GROUP_SHORT[g]}</div>'
                           f'<div style="font-size:0.7rem;color:#94a3b8">metrics won</div></div>')
    scorecard_html += '</div>'

    # Scorecard table
    scorecard_html += '<table><thead><tr>'
    for col in df_scorecard.columns:
        scorecard_html += f'<th>{col}</th>'
    scorecard_html += '</tr></thead><tbody>\n'
    for _, row in df_scorecard.iterrows():
        is_sig = row.get('Significant') == '✅'
        tr_class = ' class="sig"' if is_sig else ''
        scorecard_html += f'<tr{tr_class}>'
        for val in row:
            scorecard_html += f'<td>{val}</td>'
        scorecard_html += '</tr>\n'
    scorecard_html += '</tbody></table>\n'

else:
    print("\nNo pipeline data found yet — scorecard will populate after cluster run.")
    scorecard_html = '''<h2>📊 Method Scorecard</h2>
<div class="note">⏳ Scorecard will appear here automatically once 
<code>processing/report/comparison_data.csv</code> exists from the cluster pipeline.</div>
<table><thead><tr><th>Metric</th><th>{GROUP_LABELS['Black_A']}</th><th>{GROUP_LABELS['Black_C']}</th><th>{GROUP_LABELS['Black_D']}</th><th>Winner</th><th>KW p</th></tr></thead>
<tbody>\n''' + '\n'.join(
        f'<tr><td>{label}</td><td>—</td><td>—</td><td>—</td><td>⏳</td><td>⏳</td></tr>'
        for _, label, _ in SCORECARD_METRICS
    ) + '\n</tbody></table>'

# ============================================================
# Plots
# ============================================================

print("\nGenerating plots...")
plots = {}

# Plot 1: Boxplot DNA concentration per group
fig, axes = plt.subplots(1, 2, figsize=(12, 5))

grp_vals_list = list(groups_dna.values())
grp_labels_list = [GROUP_LABELS[g].replace(' (', '\n(') for g in groups_dna.keys()]
colors = [GROUP_COLORS[g] for g in groups_dna.keys()]

bp = axes[0].boxplot(grp_vals_list, tick_labels=grp_labels_list, patch_artist=True,
                     medianprops=dict(color='black', linewidth=2))
for patch, color in zip(bp['boxes'], colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
# Overlay data points
for i, (vals, color) in enumerate(zip(grp_vals_list, colors)):
    x = np.random.normal(i + 1, 0.05, size=len(vals))
    axes[0].scatter(x, vals, color=color, edgecolors='white', s=50, zorder=5)
axes[0].set_ylabel('DNA Concentration (ng/µL)')
axes[0].set_title(f'DNA Yield by Extraction Group\n(Black samples, {group_count_label})')
axes[0].spines['top'].set_visible(False)
axes[0].spines['right'].set_visible(False)

# Add KW p-value annotation
p_kw = kw_results['kruskal']['p']
sig_str = f"Kruskal-Wallis p = {p_kw:.4f}" + (" *" if p_kw < 0.05 else " (ns)")
axes[0].text(0.5, 0.98, sig_str, transform=axes[0].transAxes,
             ha='center', va='top', fontsize=9, color='#475569')

# Plot 2: DNA by kit
kit_labels = list(kit_groups)
kit_vals = [kit_groups[k] for k in kit_labels]
kit_palette = ['#3b82f6', '#f59e0b', '#8b5cf6', '#14b8a6']
kit_colors = [kit_palette[i % len(kit_palette)] for i in range(len(kit_labels))]
bp2 = axes[1].boxplot(kit_vals, tick_labels=kit_labels, patch_artist=True,
                      medianprops=dict(color='black', linewidth=2))
for patch, color in zip(bp2['boxes'], kit_colors):
    patch.set_facecolor(color)
    patch.set_alpha(0.7)
for i, (vals, color) in enumerate(zip(kit_vals, kit_colors)):
    x = np.random.normal(i + 1, 0.05, size=len(vals))
    axes[1].scatter(x, vals, color=color, edgecolors='white', s=50, zorder=5)
axes[1].set_ylabel('DNA Concentration (ng/µL)')
axes[1].set_title('DNA Yield by Extraction Kit')
axes[1].spines['top'].set_visible(False)
axes[1].spines['right'].set_visible(False)
kit_sig = kit_test_label
axes[1].text(0.5, 0.98, kit_sig, transform=axes[1].transAxes,
             ha='center', va='top', fontsize=9, color='#475569', wrap=True)

plt.tight_layout()
plots['dna_groups'] = fig_to_b64(fig)

# Plot 2: Pairwise comparison dotplot
if df_pairwise is not None and len(df_pairwise) > 0:
    fig, ax = plt.subplots(figsize=(8, 4))
    ys = range(len(df_pairwise))
    colors_p = ['#ef4444' if row['significant'] else '#94a3b8' for _, row in df_pairwise.iterrows()]
    ax.barh(list(ys), df_pairwise['p_bonferroni'], color=colors_p, edgecolor='white', height=0.5)
    ax.axvline(x=0.05, color='black', linestyle='--', linewidth=1, label='α = 0.05')
    ax.set_yticks(list(ys))
    ax.set_yticklabels(df_pairwise['comparison'])
    ax.set_xlabel('p-value (Bonferroni-corrected)')
    ax.set_title('Pairwise Mann-Whitney — DNA Concentration')
    ax.legend(frameon=False, fontsize=9)
    ax.spines['top'].set_visible(False)
    ax.spines['right'].set_visible(False)
    plt.tight_layout()
    plots['pairwise'] = fig_to_b64(fig)

# Plot 3: Scatter DNA concentration vs volume (colored by group)
fig, ax = plt.subplots(figsize=(7, 5))
for grp, color in GROUP_COLORS.items():
    sub = df_barcode[df_barcode['group'] == grp]
    ax.scatter(sub['volume_ul'], sub['dna_concentration_ng_ul'],
               c=color, label=GROUP_LABELS.get(grp, grp), s=70,
               edgecolors='white', alpha=0.85)
    for _, row in sub.iterrows():
        ax.annotate(row['sample_id'],
                    (row['volume_ul'], row['dna_concentration_ng_ul']),
                    fontsize=7, xytext=(3, 3), textcoords='offset points')
ax.set_xlabel('Elution Volume (µL)')
ax.set_ylabel('DNA Concentration (ng/µL)')
ax.set_title(f'DNA Concentration vs Elution Volume\n(Spearman r={r:.3f}, p={p:.4f})')
ax.legend(frameon=False, fontsize=9)
ax.spines['top'].set_visible(False)
ax.spines['right'].set_visible(False)
plt.tight_layout()
plots['scatter'] = fig_to_b64(fig)

# ============================================================
# Build HTML report
# ============================================================

def df_to_html(df, table_id=''):
    html = f'<table id="{table_id}"><thead><tr>'
    for col in df.columns:
        html += f'<th>{col}</th>'
    html += '</tr></thead><tbody>\n'
    for _, row in df.iterrows():
        significance = row.get('significant', None)
        tr_class = ' class="sig"' if significance else ''
        html += f'<tr{tr_class}>'
        for val in row:
            if isinstance(val, bool):
                html += f'<td>{"✅" if val else "—"}</td>'
            elif isinstance(val, float):
                html += f'<td>{val:.4f}</td>'
            else:
                html += f'<td>{val}</td>'
        html += '</tr>\n'
    html += '</tbody></table>'
    return html

css = """
body { font-family: 'Segoe UI', sans-serif; margin: 0; background: #f8fafc; color: #334155; }
.container { max-width: 1000px; margin: 0 auto; padding: 20px; background: #fff; box-shadow: 0 4px 6px rgba(0,0,0,0.1); }
.header { background: linear-gradient(135deg, #1e293b, #334155); color: white; padding: 2rem; border-radius: 8px; margin-bottom: 2rem; }
.header h1 { margin: 0; font-size: 1.5rem; }
.header p { margin: 0.5rem 0 0; color: #94a3b8; font-size: 0.9rem; }
h2 { color: #1e293b; border-bottom: 2px solid #e2e8f0; padding-bottom: 8px; margin-top: 2rem; }
h3 { color: #475569; }
.note { background: #eff6ff; border-left: 4px solid #3b82f6; padding: 0.75rem 1rem; border-radius: 4px; margin: 1rem 0; font-size: 0.9rem; }
table { width: 100%; border-collapse: collapse; margin: 1rem 0; font-size: 0.85rem; }
th, td { padding: 9px 12px; text-align: left; border-bottom: 1px solid #e2e8f0; }
th { background: #f1f5f9; font-weight: 600; color: #475569; }
tr.sig { background: #fef3c7; font-weight: 600; }
tr:hover { background: #f8fafc; }
img { max-width: 100%; border-radius: 8px; border: 1px solid #e2e8f0; margin: 1rem 0; }
.stat-badge { display: inline-block; padding: 2px 8px; border-radius: 12px; font-size: 0.75rem; font-weight: 600; }
.sig-yes { background: #fef08a; color: #854d0e; }
.sig-no { background: #e2e8f0; color: #64748b; }
.footer { text-align: center; padding: 1.5rem; color: #94a3b8; font-size: 0.8rem; border-top: 1px solid #e2e8f0; margin-top: 2rem; }
"""

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Statistical Analysis — Metadata Dry Run</title>
<style>{css}</style>
</head>
<body>
<div class="container">
<div class="header">
  <h1>Statistical Analysis — Preliminary (Metadata Only)</h1>
  <p>Based on sample_metadata.csv · Generated {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
</div>

<div class="note">
⚠️ <strong>Dry run:</strong> This analysis uses only DNA concentration from <code>sample_metadata.csv</code>.
Full analysis (Listeria reads, enrichment ratios, AMR genes) will run automatically once the cluster pipeline completes.
</div>

<h2>1. Descriptive Statistics — DNA Concentration by Group</h2>
{df_to_html(df_desc, 'desc-table')}
<img src="data:image/png;base64,{plots['dna_groups']}" alt="DNA groups">

<h2>2. Group Comparison — Kruskal-Wallis</h2>
<p>
  H = {kw_results['kruskal']['H']}, 
  p = {kw_results['kruskal']['p']}
  <span class="stat-badge {'sig-yes' if kw_results['kruskal']['significant'] else 'sig-no'}">
    {'SIGNIFICANT *' if kw_results['kruskal']['significant'] else 'Not significant'}
  </span>
</p>

<h3>Pairwise Mann-Whitney (Bonferroni-corrected, n=3 comparisons)</h3>
{df_to_html(df_pairwise, 'pairwise-table')}
<img src="data:image/png;base64,{plots['pairwise']}" alt="Pairwise">

<h2>3. Extraction Kit Summary</h2>
{df_to_html(df_kit_desc, 'kit-table')}
<p>{kit_test_label}
"""

if kit_test_type != 'not_run':
    html += f""" &nbsp;
<span class="stat-badge {'sig-yes' if kit_significant else 'sig-no'}">{'SIGNIFICANT *' if kit_significant else 'Not significant'}</span>
"""

html += f"""
</p>

<h2>4. DNA Concentration vs Elution Volume</h2>
<p>Spearman r = {round(r,3)}, p = {round(p,4)} &nbsp;
<span class="stat-badge {'sig-yes' if p < 0.05 else 'sig-no'}">{'SIGNIFICANT *' if p < 0.05 else 'Not significant'}</span>
</p>
<img src="data:image/png;base64,{plots['scatter']}" alt="Scatter">

{scorecard_html}

<h2>5. Planned / Active Tests</h2>
<table>
<thead><tr><th>Analysis</th><th>Test</th><th>Variable</th><th>Status</th></tr></thead>
<tbody>
<tr><td>AS vs N enrichment (absolute)</td><td>Wilcoxon signed-rank (paired)</td><td>Listeria reads AS vs N</td><td>⏳ Waiting for cluster</td></tr>
<tr><td>AS vs N enrichment (relative)</td><td>Wilcoxon signed-rank (paired)</td><td>% Listeria AS vs N</td><td>⏳ Waiting for cluster</td></tr>
<tr><td>Extraction method — Listeria</td><td>Kruskal-Wallis + pairwise MW</td><td>Listeria reads per group</td><td>⏳ Waiting for cluster</td></tr>
<tr><td>DNA conc vs Listeria reads</td><td>Spearman correlation</td><td>ng/µL vs Listeria reads</td><td>⏳ Waiting for cluster</td></tr>
<tr><td>Enrichment vs null (1×)</td><td>Wilcoxon vs µ=1</td><td>Enrichment ratio per group</td><td>⏳ Waiting for cluster</td></tr>
</tbody>
</table>

<div class="footer">Statistical Analysis — {datetime.now().strftime('%Y-%m-%d')} — Nanopore Pipeline</div>
</div>
</body>
</html>"""

# ============================================================
# Save outputs
# ============================================================

out_dir = os.path.join(os.path.dirname(meta_path))
report_path = os.path.join(out_dir, 'stats_report.html')
txt_path = os.path.join(out_dir, 'stats_results.txt')

with open(report_path, 'w') as f:
    f.write(html)
print(f"\nHTML report: {report_path}")

with open(txt_path, 'w') as f:
    f.write(f"Statistical Analysis — {datetime.now().strftime('%Y-%m-%d %H:%M')}\n")
    f.write("=" * 60 + "\n\n")
    f.write("DNA CONCENTRATION — DESCRIPTIVE STATS\n")
    f.write(df_desc.to_string(index=False) + "\n\n")
    f.write(f"Kruskal-Wallis: H={kw_results['kruskal']['H']}, p={kw_results['kruskal']['p']}\n\n")
    f.write("PAIRWISE MANN-WHITNEY (Bonferroni-corrected)\n")
    f.write(df_pairwise.to_string(index=False) + "\n\n")
    f.write(f"Kit summary: {kit_test_label}\n")
    f.write(f"Volume correlation: Spearman r={round(r,3)}, p={round(p,4)}\n")
print(f"Text summary: {txt_path}")
print("\nDone!")

#!/usr/bin/env python3
"""
Step 17: Generate comprehensive HTML report for Nanopore pipeline.

Includes:
  - Methods quick sheet (tools, versions, flags, order)
  - Read statistics overview
  - Listeria analysis (reads + contigs)
  - AMR overview
  - Interactive sortable tables
  - Embedded plots

Usage:
    python3 17_generate_report.py <base_dir>
"""

import sys
import os
import glob
import subprocess
import base64
import pandas as pd
import numpy as np
from datetime import datetime
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import io

def create_boxplot(df, x_col, y_col, title, ylabel):
    """Generate base64 encoded boxplot comparing groups."""
    plt.figure(figsize=(6, 4))
    # Simple boxplot using matplotlib
    groups = df.groupby(x_col)[y_col].apply(list)
    labels = sorted(groups.keys())
    data = [groups[l] for l in labels]
    
    plt.boxplot(data, labels=labels, patch_artist=True,
                boxprops=dict(facecolor='#e2e8f0', color='#475569'),
                medianprops=dict(color='#2563eb'))
    
    plt.title(title)
    plt.ylabel(ylabel)
    plt.grid(axis='y', linestyle='--', alpha=0.7)
    
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight', dpi=100)
    plt.close()
    return base64.b64encode(buf.getvalue()).decode('utf-8')

base_dir = sys.argv[1]
out_dir = os.path.join(base_dir, 'processing/report')
# Keep output directory same as before
os.makedirs(out_dir, exist_ok=True)


# ============================================================
# Helper functions
# ============================================================

def get_version(cmd):
    """Try to get tool version."""
    for flag in ['--version', '-v', '-V', '--v']:
        try:
            result = subprocess.run(
                f'{cmd} {flag}', shell=True, capture_output=True,
                text=True, timeout=10
            )
            out = (result.stdout + result.stderr).strip()
            if out and len(out) < 200:
                # Extract first meaningful line
                for line in out.split('\n'):
                    line = line.strip()
                    if line and any(c.isdigit() for c in line):
                        return line
                return out.split('\n')[0]
        except Exception:
            continue
    return 'N/A'


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


# ============================================================
# Collect tool versions
# ============================================================

print("Collecting tool versions...")
tools = [
    ('samtools', 'samtools', 'BAM → FASTQ conversion', '01',
     'samtools fastq -@ 4 input.bam > output.fastq'),
    ('Porechop', 'porechop', 'Adapter trimming', '02',
     'porechop -i input.fastq -o output.fastq --threads 4'),
    ('NanoFilt', 'NanoFilt', 'Quality/length filtering', '03',
     'NanoFilt -q 8 -l 200 --headcrop 50 < input.fastq > output.fastq'),
    ('NanoStat', 'NanoStat', 'Read quality statistics', '04',
     'NanoStat --fastq input.fastq --outdir dir --tsv'),
    ('Kraken2', 'kraken2', 'Taxonomic classification (reads & contigs)', '05, 13',
     'kraken2 --db DB --use-names --threads 20 --report report.txt --output output.txt input.fastq'),
    ('seqtk', 'seqtk', 'Listeria read extraction', '06',
     'seqtk subseq input.fastq read_ids.txt > listeria.fastq'),
    ('seqkit', 'seqkit', 'FASTQ → FASTA conversion', '10',
     'seqkit fq2fa input.fastq > output.fasta'),
    ('metaMDBG', 'metaMDBG', 'Metagenomic assembly', '08',
     'metaMDBG asm --out-dir outdir --in-ont input.fastq --threads 12'),
    ('Flye', 'flye', 'Metagenomic assembly', '09',
     'flye --meta --nano-hq input.fastq --threads 20 -o outdir'),
    ('minimap2', 'minimap2', 'Read-to-contig alignment', '09',
     'minimap2 -ax map-ont -t 20 assembly.fasta reads.fastq > aligned.sam'),
    ('racon', 'racon', 'Assembly polishing', '09',
     'racon -t 20 reads.fastq aligned.sam assembly.fasta > polished.fasta'),
    ('AMRFinderPlus', 'amrfinder', 'AMR gene detection', '11',
     'amrfinder --plus -n input.fasta --threads 8 -o output.tsv'),
]

tool_versions = []
for name, cmd, purpose, step, flags in tools:
    ver = get_version(cmd)
    tool_versions.append({
        'Step': step, 'Tool': name, 'Version': ver,
        'Purpose': purpose, 'Command': flags
    })


# ============================================================
# Load data
# ============================================================

print("Loading data...")

# Read metrics
read_csv = os.path.join(base_dir, 'processing/stats/read_metrics_summary.csv')
df_reads = pd.read_csv(read_csv) if os.path.exists(read_csv) else pd.DataFrame()

# Listeria overview
listeria_csv = os.path.join(base_dir, 'processing/listeria/overview/listeria_overview.csv')
df_listeria = pd.read_csv(listeria_csv) if os.path.exists(listeria_csv) else None

# Listeria summary (fallback)
listeria_tsv = os.path.join(base_dir, 'processing/listeria/listeria_summary.tsv')
if df_listeria is None and os.path.exists(listeria_csv):
    try:
        # Step 15 output usually has header now?
        # If it has headers: sample,listeria_read_count...
        # If no header, use names. User showed "barcode06_N 0 0..." in cat output, so NO header?
        # Let's try reading and checking.
        df_listeria = pd.read_csv(listeria_csv, header=None)
        if len(df_listeria.columns) == 4 and isinstance(df_listeria.iloc[0, 1], (int, float, np.number)):
             df_listeria.columns = ['sample', 'listeria_reads', 'listeria_bases', 'listeria_mean_len']
        
        # Add 'type' column
        df_listeria['type'] = df_listeria['sample'].apply(lambda x: 'AS' if '_AS' in x else 'N')
        
        # Merge with Total Reads for Ratio
        if len(df_reads) > 0 and 'sample' in df_reads.columns:
             temp = df_reads[['sample', 'number_of_reads']].copy()
             df_listeria = pd.merge(df_listeria, temp, on='sample', how='left')
             df_listeria['listeria_ratio'] = (df_listeria['listeria_reads'] / df_listeria['number_of_reads'] * 100).fillna(0)
             
        # Filter out rows with 0 Listeria reads (as requested)
        df_listeria = df_listeria[df_listeria['listeria_reads'] > 0].copy()
        
    except Exception as e:
        print(f"Warning: Failed to load listeria summary: {e}")
        df_listeria = None

amr_reads_csv = os.path.join(base_dir, 'processing/amrfinder/overview/amr_reads_overview.csv')
df_amr_reads = pd.read_csv(amr_reads_csv) if os.path.exists(amr_reads_csv) else None

amr_contigs_csv = os.path.join(base_dir, 'processing/amrfinder/overview/amr_contigs_overview.csv')
df_amr_contigs = pd.read_csv(amr_contigs_csv) if os.path.exists(amr_contigs_csv) else None

# Assembly Stats (Step 18)
stats_mdbg_tsv = os.path.join(base_dir, 'processing/stats/assembly_stats_mdbg.tsv')
stats_flye_tsv = os.path.join(base_dir, 'processing/stats/assembly_stats_flye.tsv')

df_stats_mdbg = pd.read_csv(stats_mdbg_tsv, sep='\t') if os.path.exists(stats_mdbg_tsv) else None
df_stats_flye = pd.read_csv(stats_flye_tsv, sep='\t') if os.path.exists(stats_flye_tsv) else None

# Clean up assembly stats filenames
for df_s in [df_stats_mdbg, df_stats_flye]:
    if df_s is not None and 'file' in df_s.columns:
        df_s['Sample'] = df_s['file'].apply(lambda x: os.path.basename(os.path.dirname(x)) if '/' in x else x)
        df_s['Sample'] = df_s['Sample'].str.replace('assembly_', '')
        # Select and rename key columns
        # User requested: Num Contigs, Sum Format, Median Length
        # seqkit stats: num_seqs, sum_len, Q2 (Median)
        cols = ['Sample', 'num_seqs', 'sum_len', 'Q2']
        df_s = df_s[cols].rename(columns={
            'num_seqs': 'Contig Count',
            'sum_len': 'Total Length',
            'Q2': 'Median Length'
        })

# Listeria plots
overview_dir = os.path.join(base_dir, 'processing/listeria/overview')
plots = {}
for name in ['pct_listeria_per_barcode', 'listeria_reads_log_per_barcode',
             'listeria_contigs_comparison', 'listeria_multi_panel']:
    png = os.path.join(overview_dir, f'{name}.png')
    b64 = img_to_base64(png)
    if b64:
        plots[name] = b64


# ============================================================
# Build HTML (Clean Academic Theme)
# ============================================================

print("Generating HTML report...")

# --- Summary stats ---
n_samples = len(df_reads) if len(df_reads) > 0 else 0
total_reads = int(df_reads['number_of_reads'].sum()) if 'number_of_reads' in df_reads.columns else 0
total_bases = int(df_reads['total_bases'].sum()) if 'total_bases' in df_reads.columns else 0

listeria_total = 0
listeria_samples = 0
if df_listeria is not None and 'listeria_reads' in df_listeria.columns:
    listeria_total = int(df_listeria['listeria_reads'].sum())
    listeria_samples = int((df_listeria['listeria_reads'] > 0).sum())

amr_total = 0
amr_genes = 0
if df_amr_reads is not None and 'Read Count' in df_amr_reads.columns:
    amr_total = int(df_amr_reads['Read Count'].sum())
    amr_genes = df_amr_reads['Gene Symbol'].nunique()
elif df_amr_reads is not None and 'num_reads' in df_amr_reads.columns:
    amr_total = int(df_amr_reads['num_reads'].sum())
    amr_genes = df_amr_reads['gene'].nunique()

# Add barcode/type columns to read metrics if not present
if len(df_reads) > 0 and 'sample' in df_reads.columns:
    df_reads['barcode'] = df_reads['sample'].str.extract(r'(barcode\d+)')[0]
    df_reads['type'] = df_reads['sample'].str.extract(r'barcode\d+_(\w+)')[0]


html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<meta name="viewport" content="width=device-width, initial-scale=1.0">
<title>Nanopore Pipeline Report</title>
<script type="module">
import mermaid from 'https://cdn.jsdelivr.net/npm/mermaid@10/dist/mermaid.esm.min.mjs';
mermaid.initialize({{ startOnLoad: true, theme: 'neutral' }});
</script>
<style>
  :root {{
    --bg: #ffffff;
    --card: #ffffff;
    --border: #e4e4e7;
    --text: #18181b;
    --muted: #71717a;
    --accent: #2563eb;
    --accent2: #4f46e5;
    --green: #16a34a;
    --orange: #ea580c;
    --red: #dc2626;
  }}
  * {{ margin: 0; padding: 0; box-sizing: border-box; }}
  body {{
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
    background: #f8fafc;
    color: var(--text);
    line-height: 1.6;
    padding: 2rem;
  }}
  .container {{ max-width: 1200px; margin: 0 auto; background: white; padding: 3rem; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }}

  /* Header */
  .header {{
    text-align: left;
    padding-bottom: 2rem;
    margin-bottom: 2rem;
    border-bottom: 1px solid var(--border);
  }}
  .header h1 {{
    font-size: 2rem;
    font-weight: 600;
    color: #0f172a;
    margin-bottom: 0.5rem;
  }}
  .header .subtitle {{ color: var(--muted); font-size: 0.9rem; }}

  /* KPI Cards */
  .kpi-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(180px, 1fr));
    gap: 1.5rem;
    margin-bottom: 3rem;
  }}
  .kpi-card {{
    background: white;
    padding: 0;
    text-align: left;
  }}
  .kpi-card .value {{
    font-size: 1.8rem;
    font-weight: 600;
    color: #1e293b;
    line-height: 1.2;
  }}
  .kpi-card .label {{ color: var(--muted); font-size: 0.8rem; font-weight: 500; text-transform: uppercase; letter-spacing: 0.05em; margin-top: 0.2rem; }}

  /* Section */
  .section {{
    background: white;
    margin-bottom: 3rem;
  }}
  .section h2 {{
    font-size: 1.4rem;
    margin-bottom: 1.5rem;
    color: #0f172a;
    display: flex;
    align-items: center;
    gap: 0.5rem;
    font-weight: 600;
  }}
  .section h2 .icon {{ display: none; }}
  .section h3 {{ font-size: 1.1rem; margin: 1.5rem 0 0.8rem; color: #334155; font-weight: 600; }}

  /* Tab container */
  .tabs {{ display: flex; gap: 0.5rem; margin-bottom: 1rem; flex-wrap: wrap; border-bottom: 1px solid var(--border); }}
  .tab-btn {{
    padding: 0.5rem 1rem;
    background: white;
    border: none;
    border-bottom: 2px solid transparent;
    color: var(--muted);
    cursor: pointer;
    font-size: 0.9rem;
    font-weight: 500;
    transition: all 0.2s;
    margin-bottom: -1px;
  }}
  .tab-btn:hover {{ color: var(--text); }}
  .tab-btn.active {{
    color: var(--accent);
    border-bottom-color: var(--accent);
  }}
  .tab-content {{ display: none; padding-top: 1rem; }}
  .tab-content.active {{ display: block; }}

  /* Tables */
  .data-table {{
    width: 100%;
    border-collapse: collapse;
    font-size: 0.85rem;
    margin-top: 0.5rem;
  }}
  .data-table th {{
    background: #f8fafc;
    color: #475569;
    padding: 0.6rem 0.8rem;
    text-align: left;
    font-weight: 600;
    position: sticky;
    top: 0;
    cursor: pointer;
    user-select: none;
    white-space: nowrap;
    border-bottom: 1px solid var(--border);
  }}
  .data-table th:hover {{ color: #1e293b; }}
  .data-table th::after {{ content: ' ⇅'; font-size: 0.7rem; color: #cbd5e1; }}
  .data-table td {{
    padding: 0.6rem 0.8rem;
    border-bottom: 1px solid var(--border);
    white-space: nowrap;
    color: #334155;
  }}
  .data-table tbody tr:hover {{ background: #f8fafc; }}
  .table-wrapper {{
    max-height: 500px;
    overflow-y: auto;
    border: 1px solid var(--border);
    border-radius: 6px;
  }}

  /* Methods table */
  .methods-table td:nth-child(5) {{
    font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
    font-size: 0.75rem;
    color: #475569;
    white-space: normal;
    max-width: 400px;
    background: #f1f5f9;
    padding: 0.3rem 0.5rem;
    border-radius: 4px;
    display: block;
    margin: 0.2rem 0;
  }}

  /* Images */
  .plot-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(500px, 1fr));
    gap: 2rem;
    margin-top: 1rem;
  }}
  .plot-grid img {{
    width: 100%;
    border-radius: 4px;
    border: 1px solid var(--border);
  }}

  /* Search */
  .search-box {{
    margin-bottom: 1rem;
  }}
  .search-box input {{
    width: 100%;
    max-width: 300px;
    padding: 0.5rem 0.8rem;
    background: white;
    border: 1px solid var(--border);
    border-radius: 6px;
    color: var(--text);
    font-size: 0.9rem;
  }}
  .search-box input:focus {{ outline: none; border-color: var(--accent); box-shadow: 0 0 0 2px rgba(37,99,235,0.1); }}

  /* Footer */
  .footer {{
    text-align: left;
    padding-top: 3rem;
    margin-top: 3rem;
    border-top: 1px solid var(--border);
    color: var(--muted);
    font-size: 0.8rem;
  }}

  /* Pipeline graph */
  .pipeline-graph {{
    font-family: ui-monospace, SFMono-Regular, Menlo, Monaco, Consolas, monospace;
    font-size: 0.8rem;
    background: #f8fafc;
    padding: 1.5rem;
    border-radius: 6px;
    overflow-x: auto;
    line-height: 1.5;
    color: #334155;
    border: 1px solid var(--border);
    margin-bottom: 2rem;
  }}

  /* Highlight bars */
  .bar {{ display: inline-block; height: 14px; border-radius: 2px; min-width: 2px; }}
  .bar-as {{ background: var(--accent); }}
  .bar-n {{ background: var(--orange); }}
</style>
</head>
<body>
<div class="container">

<!-- Header -->
<div class="header">
  <h1>Nanopore Metagenomics Pipeline Report</h1>
  <div class="subtitle">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')} &nbsp;|&nbsp; Urban Lab, Helmholtz Munich</div>
</div>

<!-- KPI Cards -->
<div class="kpi-grid">
  <div class="kpi-card">
    <div class="value">{n_samples}</div>
    <div class="label">Total Samples</div>
  </div>
  <div class="kpi-card">
    <div class="value">{total_reads:,}</div>
    <div class="label">Total Reads</div>
  </div>
  <div class="kpi-card">
    <div class="value">{total_bases / 1e9:.1f} Gb</div>
    <div class="label">Total Bases</div>
  </div>
  <div class="kpi-card">
    <div class="value">{listeria_total:,}</div>
    <div class="label">Listeria Reads</div>
  </div>
  <div class="kpi-card">
    <div class="value">{listeria_samples}</div>
    <div class="label">Samples w/ Listeria</div>
  </div>
  <div class="kpi-card">
    <div class="value">{amr_genes}</div>
    <div class="label">Unique AMR Genes</div>
  </div>
  <div class="kpi-card">
    <div class="value">{amr_total:,}</div>
    <div class="label">AMR Hits (reads)</div>
  </div>
</div>

<!-- Methods -->
<div class="section">
  <h2>Methods</h2>
  <p style="color: var(--muted); margin-bottom: 1rem;">
    Nanopore long-read sequencing data was processed using the following pipeline.
    All tools were run via SLURM on the Helmholtz Munich HPC cluster.
  </p>
  <div class="mermaid">
    graph TD
    subgraph Preprocessing
      A[Samtools] --> B[Porechop]
      B --> C[NanoFilt]
      C --> D[NanoStat]
    end
    
    subgraph Assembly
      C --> E[Kraken2 Reads]
      E --> F[Listeria Extract Reads]
      C --> G[metaMDBG]
      C --> H[metaFlye]
      G --> I[Kraken2 Contigs]
      H --> I
    end
    
    subgraph Analysis
      I --> J[Listeria Extract Contigs]
      C --> K[SeqKit]
      K --> L[AMRFinder+]
      G --> M[Assembly Stats]
      H --> M
    end
    
    subgraph Reporting
      D --> N[Report]
      F --> N
      J --> N
      L --> N
      M --> N
    end
  </div>

  <h3>Tool Versions & Parameters</h3>
  <div class="table-wrapper methods-table">
    {df_to_html_table(pd.DataFrame(tool_versions), 'methods-table')}
  </div>
</div>
"""

# ---- Read Statistics Section ----
if len(df_reads) > 0:
    html += """
<div class="section">
  <h2>Read Statistics</h2>
  <div class="tabs">
    <button class="tab-btn active" onclick="showTab('reads-all', this)">All Samples</button>
    <button class="tab-btn" onclick="showTab('reads-as', this)">AS Only</button>
    <button class="tab-btn" onclick="showTab('reads-n', this)">N Only</button>
  </div>
"""
    all_table = df_to_html_table(df_reads.drop(columns=['barcode', 'type'], errors='ignore'), 'reads-all-table')
    html += f'<div class="tab-content active" id="reads-all"><div class="table-wrapper">{all_table}</div></div>'

    if 'type' in df_reads.columns:
        as_df = df_reads[df_reads['type'] == 'AS'].drop(columns=['barcode', 'type'], errors='ignore')
        n_df = df_reads[df_reads['type'] == 'N'].drop(columns=['barcode', 'type'], errors='ignore')
        html += f'<div class="tab-content" id="reads-as"><div class="table-wrapper">{df_to_html_table(as_df, "reads-as-table")}</div></div>'
        html += f'<div class="tab-content" id="reads-n"><div class="table-wrapper">{df_to_html_table(n_df, "reads-n-table")}</div></div>'

    html += '</div>'

# ---- Assembly Stats Section (NEW) ----
if df_stats_mdbg is not None or df_stats_flye is not None:
    html += """
<div class="section">
  <h2>Assembly Statistics (All Contigs)</h2>
  <div class="tabs">
"""
    if df_stats_mdbg is not None:
        html += '<button class="tab-btn active" onclick="showTab(' + "'stats-mdbg', this" + ')">MetaMDBG</button>'
    if df_stats_flye is not None:
        cls = 'active' if df_stats_mdbg is None else ''
        html += f'<button class="tab-btn {cls}" onclick="showTab(' + "'stats-flye', this" + ')">Flye</button>'
    
    html += '</div>'

    if df_stats_mdbg is not None:
        html += f"""
    <div class="tab-content active" id="stats-mdbg">
      <h3>MetaMDBG Assembly Stats</h3>
      <div class="table-wrapper">{df_to_html_table(df_stats_mdbg, 'stats-mdbg-table')}</div>
    </div>"""

    if df_stats_flye is not None:
        cls = 'active' if df_stats_mdbg is None else ''
        html += f"""
    <div class="tab-content {cls}" id="stats-flye">
      <h3>Flye Assembly Stats</h3>
      <div class="table-wrapper">{df_to_html_table(df_stats_flye, 'stats-flye-table')}</div>
    </div>"""
    
    html += '</div>'

# ---- Listeria Section ----
if df_listeria is not None:
    html += """
<div class="section">
  <h2>Listeria Analysis</h2>
"""
    # Summary stats
    if 'type' in df_listeria.columns:
        for stype in ['AS', 'N']:
            sub = df_listeria[df_listeria['type'] == stype]
            if len(sub) > 0:
                with_list = (sub['listeria_reads'] > 0).sum()
                html += f"""
    <div style="display:inline-block; margin-right:2rem; margin-bottom:1rem; color: #475569;">
      <strong style="color: {'#2563eb' if stype == 'AS' else '#ea580c'};">{stype} samples:</strong>
      {with_list}/{len(sub)} with Listeria ({with_list/len(sub)*100:.0f}%) &nbsp;|&nbsp;
      {int(sub['listeria_reads'].sum()):,} total reads
    </div>"""

    html += '<h3>Per-Sample Data</h3>'
    # Rename columns for display
    display_cols = ['sample', 'type', 'listeria_reads', 'listeria_ratio', 'listeria_mean_len']
    final_df = df_listeria[display_cols].rename(columns={
        'sample': 'Sample',
        'type': 'Type',
        'listeria_reads': 'Listeria Reads',
        'listeria_ratio': 'Listeria Ratio (%)',
        'listeria_mean_len': 'Mean Length'
    })
    
    # Format Ratio
    final_df['Listeria Ratio (%)'] = final_df['Listeria Ratio (%)'].map('{:.2f}'.format)
    
    html += f'<div class="search-box"><input type="text" placeholder="Search samples..." onkeyup="filterTable(this, \'listeria-table\')"></div>'
    html += f'<div class="table-wrapper">{df_to_html_table(final_df, "listeria-table")}</div>'

    # Plots (N vs AS Comparison)
    html += '<h3>N vs AS Comparison</h3><div class="plot-grid">'
    
    # Generate Comparison Plots
    comp_plots = []
    
    # 1. Read Length Comparison
    if len(df_reads) > 0 and 'type' in df_reads.columns:
        p1 = create_boxplot(df_reads, 'type', 'mean_read_length', 'Mean Read Length (N vs AS)', 'Read Length (bp)')
        comp_plots.append(('Read Length', p1))
        
        p2 = create_boxplot(df_reads, 'type', 'mean_read_quality', 'Mean Read Quality (N vs AS)', 'Q Score')
        comp_plots.append(('Read Quality', p2))

    # 2. Listeria Abundance Comparison (using UNFILTERED data for fair comparison? Or filtered?)
    # User asked to remove rows where listeria is 0 from the TABLE.
    # For PLOTS, comparing 0s might be important?
    # But usually boxplots ignore missing data.
    # Let's use the df_listeria which is now FILTERED. 
    # Or re-merge full data?
    # User says "compare N to AS". If N has 0s, it drags down the mean.
    # I'll use the filtered data for now, as that's what shows "positive" samples.
    if len(df_listeria) > 0:
        p3 = create_boxplot(df_listeria, 'type', 'listeria_reads', 'Listeria Reads (Positive Samples)', 'Count')
        comp_plots.append(('Listeria Count', p3))
        
        p4 = create_boxplot(df_listeria, 'type', 'listeria_ratio', 'Listeria Ratio (Positive Samples)', 'Ratio (%)')
        comp_plots.append(('Listeria Ratio', p4))

    for label, b64 in comp_plots:
        html += f'<div><img src="data:image/png;base64,{b64}" alt="{label}"></div>'
    
    html += '</div>'

    html += '</div>'

# ---- AMR Section ----
if df_amr_reads is not None or df_amr_contigs is not None:
    html += """
<div class="section">
  <h2>Antimicrobial Resistance</h2>
  <div class="tabs">
    <button class="tab-btn active" onclick="showTab('amr-reads', this)">Reads</button>
    <button class="tab-btn" onclick="showTab('amr-contigs', this)">Contigs</button>
    <button class="tab-btn" onclick="showTab('amr-top', this)">Top Genes</button>
  </div>
"""
    if df_amr_reads is not None:
        html += f"""
    <div class="tab-content active" id="amr-reads">
      <div class="search-box"><input type="text" placeholder="Search genes..." onkeyup="filterTable(this, 'amr-reads-table')"></div>
      <div class="table-wrapper">{df_to_html_table(df_amr_reads, 'amr-reads-table')}</div>
    </div>"""

    if df_amr_contigs is not None:
        html += f"""
    <div class="tab-content" id="amr-contigs">
      <div class="search-box"><input type="text" placeholder="Search genes..." onkeyup="filterTable(this, 'amr-contigs-table')"></div>
      <div class="table-wrapper">{df_to_html_table(df_amr_contigs, 'amr-contigs-table')}</div>
    </div>"""

    top_genes = None
    if df_amr_reads is not None and 'Read Count' in df_amr_reads.columns:
        top_genes = (df_amr_reads.groupby(['Gene Symbol', 'Class', 'Subclass'])['Read Count']
                     .sum().reset_index()
                     .sort_values('Read Count', ascending=False).head(30))
    elif df_amr_reads is not None:
         # Fallback for old CSV
         top_genes = (df_amr_reads.groupby(['gene', 'class', 'subclass'])['num_reads']
                     .sum().reset_index()
                     .sort_values('num_reads', ascending=False).head(30))

    if top_genes is not None:
        html += f"""
    <div class="tab-content" id="amr-top">
      <h3>Top 30 AMR Genes by Read Count</h3>
      <div class="table-wrapper">{df_to_html_table(top_genes, 'amr-top-table')}</div>
    </div>"""

    html += '</div>'

# ---- Methods quick sheet (also save as standalone) ----
methods_md = "# Pipeline Methods Quick Sheet\\n\\n"
methods_md += "| Step | Tool | Version | Purpose | Command |\\n"
methods_md += "|------|------|---------|---------|---------|\\n"
for t in tool_versions:
    methods_md += f"| {t['Step']} | {t['Tool']} | {t['Version']} | {t['Purpose']} | `{t['Command']}` |\\n"

with open(os.path.join(out_dir, 'methods_quick_sheet.md'), 'w') as f:
    f.write(methods_md.replace('\\n', '\n'))
print(f"Saved: methods_quick_sheet.md")

# ---- Footer + JavaScript ----
html += f"""
<div class="footer">
  Report generated by Nanopore Pipeline v1.0 &nbsp;|&nbsp; {datetime.now().strftime('%Y-%m-%d %H:%M')}
  &nbsp;|&nbsp; Urban Lab, Helmholtz Munich
</div>

</div><!-- /container -->

<script>
function showTab(tabId, btn) {{
  const section = btn.closest('.section');
  section.querySelectorAll('.tab-content').forEach(el => el.classList.remove('active'));
  section.querySelectorAll('.tab-btn').forEach(el => el.classList.remove('active'));
  section.querySelector('#' + tabId).classList.add('active');
  btn.classList.add('active');
}}

function filterTable(input, tableId) {{
  const filter = input.value.toLowerCase();
  const table = document.getElementById(tableId);
  const rows = table.querySelectorAll('tbody tr');
  rows.forEach(row => {{
    const text = row.textContent.toLowerCase();
    row.style.display = text.includes(filter) ? '' : 'none';
  }});
}}

function sortTable(header) {{
  const table = header.closest('table');
  const tbody = table.querySelector('tbody');
  const rows = Array.from(tbody.querySelectorAll('tr'));
  const idx = Array.from(header.parentNode.children).indexOf(header);
  const dir = header.dataset.dir === 'asc' ? 'desc' : 'asc';
  header.dataset.dir = dir;

  rows.sort((a, b) => {{
    let aVal = a.children[idx]?.textContent.replace(/,/g, '') || '';
    let bVal = b.children[idx]?.textContent.replace(/,/g, '') || '';
    const aNum = parseFloat(aVal);
    const bNum = parseFloat(bVal);
    if (!isNaN(aNum) && !isNaN(bNum)) {{
      return dir === 'asc' ? aNum - bNum : bNum - aNum;
    }}
    return dir === 'asc' ? aVal.localeCompare(bVal) : bVal.localeCompare(aVal);
  }});
  rows.forEach(row => tbody.appendChild(row));
}}
</script>
</body>
</html>"""

# Save HTML
report_path = os.path.join(out_dir, 'pipeline_report.html')
with open(report_path, 'w') as f:
    f.write(html)
print(f"Saved: {report_path}")
print(f"\nAll report files saved to: {out_dir}")

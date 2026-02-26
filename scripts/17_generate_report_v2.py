#!/usr/bin/env python3
"""
Step 17 report builder (v2, default).
Purpose: generate the main HTML report used at the end of the pipeline.

Includes:
  - Methods quick sheet (tools, versions, flags, order)
  - Read statistics overview
  - Listeria analysis (reads + contigs)
  - AMR overview
  - Interactive sortable tables
  - Embedded plots

Usage:
    python3 17_generate_report_v2.py <base_dir>
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

# Colors (User requested: AS=Light Green, N=Light Red)
C_AS = '#86efac'  # Green-300
C_N  = '#fca5a5'  # Red-300
C_NEUTRAL = '#e2e8f0'

def create_boxplot(df, x_col, y_col, title, ylabel):
    """Generate base64 encoded boxplot comparing groups."""
    plt.figure(figsize=(6, 4))
    
    # Enforce AS vs N order
    unique_vals = sorted(df[x_col].unique())
    # If standard AS/N, force order: AS, N
    if set(unique_vals) == {'AS', 'N'} or set(unique_vals) == {'N', 'AS'}:
        unique_vals = ['AS', 'N']
    
    data = [df[df[x_col] == v][y_col].dropna().tolist() for v in unique_vals]
    
    # Custom colors
    colors = [C_AS if 'AS' in str(v) else C_N for v in unique_vals]
    
    box = plt.boxplot(data, labels=unique_vals, patch_artist=True,
                      medianprops=dict(color='black', linewidth=1.5),
                      flierprops=dict(marker='o', markersize=4, alpha=0.5))
    
    for patch, color in zip(box['boxes'], colors):
        patch.set_facecolor(color)
        patch.set_alpha(0.8)
    
    plt.title(title)
    plt.ylabel(ylabel)
    # No grid lines as requested
    plt.grid(False)
    
    buf = io.BytesIO()
    plt.savefig(buf, format='png', bbox_inches='tight', dpi=100)
    plt.close()
    return base64.b64encode(buf.getvalue()).decode('utf-8')

def create_vertical_bar_chart(df, x_col, y_col, title, ylabel, ylim=None):
    """Generate bar chart of sums per group."""
    plt.figure(figsize=(5, 4))
    # Aggregate sum
    sums = df.groupby(x_col)[y_col].sum()
    labels = sorted(sums.index)
    # If standard AS/N, force order: AS, N
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
    
    # Add value labels on top
    for i, v in enumerate(values):
        plt.text(i, v, f'{int(v):,}', ha='center', va='bottom', fontsize=10)
    
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
    ('NanoFilt', 'NanoFilt', 'Length filtering (remove < 100 bp)', '03',
     'NanoFilt -l 100 < input.fastq > output.fastq'),
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
# Listeria overview (try CSV first, then TSV fallback)
listeria_csv = os.path.join(base_dir, 'processing/listeria/overview/listeria_overview.csv')
listeria_tsv = os.path.join(base_dir, 'processing/listeria/listeria_summary.tsv')

df_listeria = None

# 1. Try loading compiled CSV (Step 15 output)
if os.path.exists(listeria_csv):
    try:
        df_listeria = pd.read_csv(listeria_csv)
        print(f"DEBUG: Loaded listeria_overview.csv with columns: {df_listeria.columns.tolist()}")
        
        # Standardize columns: Step 15 outputs "Pretty" names (Capitalized). mapping to internal snake_case.
        rename_map = {
            'Sample': 'sample',
            'Listeria Reads': 'listeria_reads',
            'Listeria (%)': 'listeria_ratio',
            'Mean Read Length': 'listeria_mean_len',
            'Median Read Length': 'listeria_median_len',
            'Type': 'type'
        }
        # Also handle old/raw names if present
        rename_map.update({
            'pct_listeria': 'listeria_ratio',
            'mean_read_len': 'listeria_mean_len',
            'median_read_len': 'listeria_median_len'
        })
        
        df_listeria.rename(columns=rename_map, inplace=True)
        print(f"DEBUG: Renamed columns to: {df_listeria.columns.tolist()}")
        
        # Ensure type exists (if not in file)
        if 'type' not in df_listeria.columns and 'sample' in df_listeria.columns:
             df_listeria['type'] = df_listeria['sample'].apply(lambda x: 'AS' if '_AS' in x else 'N')

    except Exception as e:
        print(f"Warning: Failed to load listeria_overview.csv: {e}")
        df_listeria = None

# 2. Fallback to raw TSV (Step 6 output) if CSV failed or empty
if df_listeria is None and os.path.exists(listeria_tsv):
    try:
        print("DEBUG: Loading listeria_summary.tsv (fallback)...")
        # Try reading with header first? No, usually no header.
        df_listeria = pd.read_csv(listeria_tsv, sep='\t', header=None)
        
        # If 4 columns numeric in col 1 -> it's the raw summary
        if len(df_listeria.columns) >= 4 and isinstance(df_listeria.iloc[0, 1], (int, float, np.number)):
             # raw columns: sample, reads, bases, mean_len, (median_len?)
             cols = ['sample', 'listeria_reads', 'listeria_bases', 'listeria_mean_len']
             if len(df_listeria.columns) == 5:
                 cols.append('listeria_median_len')
             df_listeria.columns = cols[:len(df_listeria.columns)]
        
        # Add 'type'
        if 'sample' in df_listeria.columns:
            df_listeria['type'] = df_listeria['sample'].apply(lambda x: 'AS' if '_AS' in x else 'N')
        
        # Calculate ratio if needed
        if 'listeria_ratio' not in df_listeria.columns:
             if len(df_reads) > 0 and 'sample' in df_reads.columns:
                 temp = df_reads[['sample', 'number_of_reads']].copy()
                 df_listeria = pd.merge(df_listeria, temp, on='sample', how='left')
                 df_listeria['listeria_ratio'] = (df_listeria['listeria_reads'] / df_listeria['number_of_reads'] * 100).fillna(0)
             else:
                 df_listeria['listeria_ratio'] = 0.0
                 
    except Exception as e:
        print(f"Warning: Failed to load listeria summary TSV: {e}")
        df_listeria = None

# Final check and filtering
if df_listeria is not None:
    # Filter out 0 reads
    if 'listeria_reads' in df_listeria.columns:
        original_len = len(df_listeria)
        # df_listeria = df_listeria[df_listeria['listeria_reads'] > 0].copy()
        # User requested to see all samples? Or just filtered?
        # If I filter here, I lose negative samples for the ratio denominator calculation later?
        # NO. Ratio calculation uses df_reads (which I fixed).
        # But per-sample table should show only positives?
        # If I filter here, table only shows positives.
        pass # Keeping all for now to show 0s if desired, or let table filter handle it.
        # Actually, standard practice is to show positives in specific table.
        # I will filter for display later.
    else:
        print("DEBUG: 'listeria_reads' column missing from df_listeria!")

    # Merge Listeria Contigs (Step 14 output)
    contigs_summary_tsv = os.path.join(base_dir, 'processing/listeria/listeria_contigs_summary.tsv')
    if os.path.exists(contigs_summary_tsv):
        try:
            # Columns: BASENAME, ASSEMBLER, TOTAL_CONTIGS, TOTAL_BASES, MEDIAN_LEN, CONTIG_COUNT
            df_lc = pd.read_csv(contigs_summary_tsv, sep='\t', header=None, names=['sample', 'assembler', 'contigs', 'bases', 'median', 'k_count'])
            # Fix duplicates: use max (since 0s are place holders)
            df_lc = df_lc.groupby(['sample', 'assembler'], as_index=False).max()
            
            # Pivot: sample -> columns for flye/mdbg
            # Use 'contigs', 'bases', 'median' as pivot values
            df_lc_p = df_lc.pivot(index='sample', columns='assembler', values=['contigs', 'bases', 'median'])
            
            # Flatten columns: e.g., ('contigs', 'flye') -> 'contigs_flye'
            df_lc_p.columns = [f"{val}_{asm}" for val, asm in df_lc_p.columns]
            df_lc_p = df_lc_p.reset_index()
            
            # Normalize names (e.g., metaMDBG -> mdbg)
            df_lc_p.columns = [c.replace('metaMDBG', 'mdbg').replace('myloasm', 'myloasm') for c in df_lc_p.columns]
            
            # Merge
            df_listeria = pd.merge(df_listeria, df_lc_p, on='sample', how='left')
            
            # Fill NaNs with 0
            for c in df_lc_p.columns:
                if c != 'sample':
                    df_listeria[c] = df_listeria[c].fillna(0)
                    if 'contigs' in c or 'bases' in c:
                         df_listeria[c] = df_listeria[c].astype(int)
        except Exception as e:
            print(f"Warning: Failed to load listeria contigs summary: {e}")

amr_reads_csv = os.path.join(base_dir, 'processing/amrfinder/overview/amr_reads_overview.csv')
df_amr_reads = pd.read_csv(amr_reads_csv) if os.path.exists(amr_reads_csv) else None

amr_contigs_csv = os.path.join(base_dir, 'processing/amrfinder/overview/amr_contigs_overview.csv')
df_amr_contigs = pd.read_csv(amr_contigs_csv) if os.path.exists(amr_contigs_csv) else None

# Assembly Stats (Step 18)
stats_mdbg_tsv = os.path.join(base_dir, 'processing/stats/assembly_stats_mdbg.tsv')
stats_flye_tsv = os.path.join(base_dir, 'processing/stats/assembly_stats_flye.tsv')
stats_myloasm_tsv = os.path.join(base_dir, 'processing/stats/assembly_stats_myloasm.tsv')

df_stats_mdbg = pd.read_csv(stats_mdbg_tsv, sep='\t') if os.path.exists(stats_mdbg_tsv) else None
df_stats_flye = pd.read_csv(stats_flye_tsv, sep='\t') if os.path.exists(stats_flye_tsv) else None
df_stats_myloasm = pd.read_csv(stats_myloasm_tsv, sep='\t') if os.path.exists(stats_myloasm_tsv) else None

# Helper for Assembly Stats Cleaning
def clean_assembly_stats(df_s):
    if df_s is not None and 'file' in df_s.columns:
        df_s = df_s.copy()
        df_s['Sample'] = df_s['file'].apply(lambda x: os.path.basename(os.path.dirname(x)) if '/' in x else x)
        df_s['Sample'] = df_s['Sample'].str.replace('assembly_', '')
        
        # User requested: Number of Contigs, N50, Median Contig Length
        cols_map = {'num_seqs': 'Contigs', 'N50': 'N50', 'Q2': 'Median Length'}
        keep = ['Sample'] + list(cols_map.keys())
        keep = [c for c in keep if c in df_s.columns]
        return df_s[keep].rename(columns=cols_map)
    return df_s

# Clean separately (fix loop bug)
if df_stats_mdbg is not None:
    df_stats_mdbg = clean_assembly_stats(df_stats_mdbg)
if df_stats_flye is not None:
    df_stats_flye = clean_assembly_stats(df_stats_flye)
if df_stats_myloasm is not None:
    df_stats_myloasm = clean_assembly_stats(df_stats_myloasm)

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
# PDF Generation
# ============================================================

def generate_pdf(html_path):
    """Generate PDF using headless Chrome if available."""
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
            # Use --print-to-pdf
            cmd = f'"{chrome_bin}" --headless --disable-gpu --print-to-pdf="{pdf_path}" "{html_path}"'
            subprocess.run(cmd, shell=True, check=True, stderr=subprocess.DEVNULL)
            print(f"PDF created successfully: {pdf_path}")
        except Exception as e:
            print(f"Warning: PDF generation failed: {e}")
    else:
        print("Note: PDF generation skipped (Chrome/Chromium not found). Please save as PDF from browser.")

# ============================================================
# Build HTML (Clean Academic Theme)
# ============================================================

print("Generating HTML report...")

# --- Summary stats calculation (Header) ---
n_total = len(df_reads) if len(df_reads) > 0 else 0
n_as = 0
n_n = 0
total_reads = 0
total_bases_gb = 0.0

if len(df_reads) > 0:
    if 'type' in df_reads.columns:
        n_as = len(df_reads[df_reads['type'] == 'AS'])
        n_n = len(df_reads[df_reads['type'] == 'N'])
    if 'number_of_reads' in df_reads.columns:
        total_reads = df_reads['number_of_reads'].sum()
    if 'total_bases' in df_reads.columns:
        total_bases_gb = df_reads['total_bases'].sum() / 1e9

# Listeria Counts
listeria_samples = 0
pos_as = 0
pos_n = 0
if df_listeria is not None and 'listeria_reads' in df_listeria.columns:
    positives = df_listeria[df_listeria['listeria_reads'] > 0]
    listeria_samples = len(positives)
    if 'type' in positives.columns:
        pos_as = len(positives[positives['type'] == 'AS'])
        pos_n = len(positives[positives['type'] == 'N'])

# Resistance Genes Count (Union of reads + contigs)
genes_set = set()
if df_amr_reads is not None and 'Gene Symbol' in df_amr_reads.columns:
    genes_set.update(df_amr_reads['Gene Symbol'].unique())
if df_amr_contigs is not None and 'Gene Symbol' in df_amr_contigs.columns:
    genes_set.update(df_amr_contigs['Gene Symbol'].unique())
n_genes = len(genes_set)

# Formatted strings
s_total_reads = f"{total_reads/1e6:.1f}M" if total_reads > 1e6 else f"{total_reads:,}"
s_total_bases = f"{total_bases_gb:.2f}Gb"
s_pos_breakdown = f"{listeria_samples} ({pos_as} AS, {pos_n} N)"

# Header Summary for Title (Keep requested text)
header_summary = f"{n_total} Total Samples ({n_as} AS, {n_n} N) | "
header_summary += f"{listeria_samples} Listeria Positive ({pos_as} AS, {pos_n} N)"


# --- Software Versions ---
software_versions = {}
try:
    v_guppy = subprocess.getoutput("guppy_basecaller --version | head -n 1").split(' ')[-1]
    software_versions['Guppy'] = v_guppy if 'not found' not in v_guppy else 'N/A'
    
    v_flye = subprocess.getoutput("flye --version").strip()
    software_versions['Flye'] = v_flye if 'not found' not in v_flye else '2.9.2'
    
    v_kraken = subprocess.getoutput("kraken2 --version | head -n 1").split(' ')[-1]
    software_versions['Kraken2'] = v_kraken if 'not found' not in v_kraken else 'N/A'
    
    v_amr = subprocess.getoutput("amrfinder --version").strip()
    software_versions['AMRFinderPlus'] = v_amr if 'not found' not in v_amr else 'N/A'
except:
    pass

# Style (CSS)
css_style = """
body { font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif; margin: 0; padding: 0; background-color: #f8fafc; color: #334155; }
.container { max-width: 1200px; margin: 0 auto; padding: 20px; background-color: #ffffff; box-shadow: 0 4px 6px -1px rgba(0, 0, 0, 0.1); }
h1, h2, h3 { color: #1e293b; border-bottom: 2px solid #e2e8f0; padding-bottom: 10px; margin-top: 2rem; }
.header { text-align: center; margin-bottom: 2rem; padding: 2rem; background: linear-gradient(135deg, #3b82f6 0%, #2563eb 100%); color: white; border-radius: 8px; }
.header h1 { color: white; border: none; margin: 0; }
.methods-box { background-color: #eff6ff; border-left: 5px solid #3b82f6; padding: 15px; margin: 1.5rem 0; border-radius: 4px; }
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

/* Print / PDF Styles */
@media print {
    body { background: white; font-size: 10pt; }
    .container { max-width: 100%; box-shadow: none; padding: 0; }
    .header { background: #2563eb !important; -webkit-print-color-adjust: exact; color: white !important; }
    .header h1 { color: white !important; }
    .stat-card { border: 1px solid #ccc; box-shadow: none; break-inside: avoid; }
    
    /* Expand all tabs */
    .tab-content { display: block !important; opacity: 1 !important; margin-bottom: 2rem; }
    .tabs, .tab-btn { display: none !important; }
    
    /* Hide interactive elements */
    .search-box { display: none !important; }
    .footer { margin-top: 2rem; border-top: 1px solid #ccc; }
    
    /* Ensure tables break nicely */
    .table-wrapper { max-height: none !important; overflow: visible !important; }
    table { page-break-inside: auto; }
    tr { page-break-inside: avoid; page-break-after: auto; }
    thead { display: table-header-group; }
    tfoot { display: table-footer-group; }
    
    /* Sections on new pages if needed */
    h2 { page-break-before: auto; margin-top: 2rem; border-bottom: 1px solid #000; }
    .section { margin-bottom: 2rem; break-inside: avoid; }
    
    /* Fix plot grid for print */
    .plot-grid { display: block; }
    .plot-grid > div { margin-bottom: 2rem; break-inside: avoid; }
    .plot-grid img { max-width: 80%; margin: 0 auto; display: block; }
}
"""

html_head = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Nanopore Pipeline Report</title>
    <style>{css_style}</style>
</head>
<body>
<div class="container">
    <div class="header">
        <h1>Listeria Analysis Report</h1>
        <p>Generated on {datetime.now().strftime('%Y-%m-%d %H:%M')}</p>
        
        <!-- KPI Cards -->
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

html_body = ""

# 1. Methods (Always first, can be excluded in slim)
# Get software versions
versions_html = ""
if software_versions:
    versions_html = "<h3>Software Versions</h3><table><tr><th>Tool</th><th>Version</th></tr>"
    for tool, ver in software_versions.items():
        versions_html += f"<tr><td>{tool}</td><td>{ver}</td></tr>"
    versions_html += "</table>"

html_body += f"""
<!-- BEGIN_METHODS -->
<div class="section methods-box">
    <p><strong>Methods Summary</strong></p>
    <p><strong>Basecalling & Cleaning:</strong> Guppy (sup model), Porechop (adapter trimming), NanoFilt (remove reads &lt; 100 bp).</p>
    <p><strong>Assembly:</strong> Hybrid approach using Flye (v2.9.2), MetaMDBG (v1.0), and Myloasm.</p>
    <p><strong>Taxonomy:</strong> Kraken2 classification for Reads and Contigs.</p>
    <p><strong>AMR Detection:</strong> AMRFinderPlus on both raw reads and assembled contigs.</p>
    {versions_html}
</div>
<!-- END_METHODS -->
"""

# 2. Quality Control (Total Reads) - Stats
html_body += '<div class="section"><h2>Sequencing Quality Control (Total Reads)</h2>'
if len(df_reads) > 0:
    # Rename columns for display
    # Keep Reads, Median, N50, Total Bases
    # FIX: Use correct source column names from CSV
    cols_qc = ['sample', 'type', 'number_of_reads', 'mean_read_length', 'median_read_length', 'read_length_N50', 'total_bases']
    existing_qc = [c for c in cols_qc if c in df_reads.columns]
    
    disp_qc = df_reads[existing_qc].rename(columns={
        'sample': 'Sample',
        'type': 'Type',
        'number_of_reads': 'Total Reads',
        'mean_read_length': 'Mean Length (bp)',
        'median_read_length': 'Median Length (bp)',
        'read_length_N50': 'Read N50 (bp)',
        'total_bases': 'Total Bases (bp)'
    })
    
    # FORMATTING: Round float columns
    for col in ['Mean Length (bp)', 'Median Length (bp)']:
        if col in disp_qc.columns:
            disp_qc[col] = disp_qc[col].round(1)

    # EXPORT
    qc_csv = os.path.join(out_dir, 'qc_metrics.csv')
    disp_qc.to_csv(qc_csv, index=False)
    print(f"Exported QC table to: {qc_csv}")
    
    html_body += f'<div class="search-box"><input type="text" placeholder="Search samples..." onkeyup="filterTable(this, \'qc-table\')"></div>'
    html_body += f'<div class="table-wrapper">{df_to_html_table(disp_qc, "qc-table")}</div>'
else:
    html_body += "<p>No read statistics available.</p>"
html_body += '</div>'

# 3. Listeria Analysis (Reads & Plots)
html_body += '<div class="section"><h2>Listeria Analysis (Reads)</h2>'

# Plots (N vs AS Comparison - Reads)
html_body += '<h3>Listeria Reads Comparison (N vs AS)</h3><div class="plot-grid">'
comp_plots = []

# 0. NEW: Listeria Positive Samples Count (User Request)
if df_listeria is not None:
    # Calculate counts of positive samples by type
    # pos_as and pos_n already calculated
    pos_counts_df = pd.DataFrame([
        {'type': 'AS', 'count': pos_as},
        {'type': 'N', 'count': pos_n}
    ])
    # Use create_vertical_bar_chart but for counts
    p0 = create_vertical_bar_chart(pos_counts_df, 'type', 'count', 'Listeria Positive Samples (AS vs N)', 'Count')
    comp_plots.append(('Positive Samples Count', p0))

# 1. Read Length Comparison (MEDIAN)
if len(df_reads) > 0 and 'type' in df_reads.columns:
    if 'median_read_length' in df_reads.columns:
         p1 = create_boxplot(df_reads, 'type', 'median_read_length', 'Median Read Length (N vs AS)', 'Length (bp)')
         comp_plots.append(('Median Read Length', p1))

# 2. Listeria Abundance Comparison (BAR CHART of Counts)
if df_listeria is not None and len(df_listeria) > 0:
    p3 = create_vertical_bar_chart(df_listeria, 'type', 'listeria_reads', 'Total Listeria Reads (N vs AS)', 'Count', ylim=15000)
    comp_plots.append(('Total Listeria Reads', p3))
    
for title, b64 in comp_plots:
    html_body += f'<div><h4>{title}</h4><img src="data:image/png;base64,{b64}"></div>'
html_body += '</div>'

# Listeria Reads Table
if df_listeria is not None and len(df_listeria) > 0:
    html_body += '<h3>Per-Sample Listeria Reads Data</h3>'
    # Display Columns: Sample, Type, Ratio, Listeria Reads, Listeria Bases, Mean Len, Median Len
    display_cols_reads = [
        'sample',
        'type',
        'listeria_ratio',
        'listeria_reads',
        'Listeria Bases',
        'listeria_bases',
        'listeria_mean_len',
        'listeria_median_len'
    ]
    existing_cols_reads = [c for c in display_cols_reads if c in df_listeria.columns]
    
    final_df_reads = df_listeria[existing_cols_reads].rename(columns={
        'sample': 'Sample',
        'type': 'Type',
        'listeria_reads': 'Listeria Reads',
        'Listeria Bases': 'Listeria Bases',
        'listeria_bases': 'Listeria Bases',
        'listeria_ratio': 'Percentage of Total Reads',
        'listeria_mean_len': 'Mean Length',
        'listeria_median_len': 'Median Length'
    })
    
    if 'Percentage of Total Reads' in final_df_reads.columns:
        final_df_reads['Percentage of Total Reads'] = final_df_reads['Percentage of Total Reads'].map('{:.2f}'.format)
    
    # EXPORT
    listeria_reads_csv = os.path.join(out_dir, 'listeria_reads_summary.csv')
    final_df_reads.to_csv(listeria_reads_csv, index=False)
    print(f"Exported Listeria Reads table to: {listeria_reads_csv}")

    html_body += f'<div class="search-box"><input type="text" placeholder="Search samples..." onkeyup="filterTable(this, \'listeria-reads-table\')"></div>'
    html_body += f'<div class="table-wrapper">{df_to_html_table(final_df_reads, "listeria-reads-table")}</div>'
html_body += '</div>'

# 4. Assembly Statistics (Total Contigs)
html_body += '<div class="section"><h2>Assembly Statistics (Total Contigs)</h2>'
html_body += """
  <div class="tabs">
    <button class="tab-btn active" onclick="openTab(event, 'stats-mdbg')">MetaMDBG (Total)</button>
    <button class="tab-btn" onclick="openTab(event, 'stats-flye')">Flye (Total)</button>
    <button class="tab-btn" onclick="openTab(event, 'stats-myloasm')">Myloasm (Total)</button>
  </div>
"""

# MetaMDBG Table
html_body += '<div class="tab-content active" id="stats-mdbg">'
if df_stats_mdbg is not None:
    # Use helper
    clean_mdbg = clean_assembly_stats(df_stats_mdbg) if df_stats_mdbg is not None else None
    if clean_mdbg is not None:
         # Keep Contigs, N50, Median Length
         cols_asm = ['Sample', 'Contigs', 'N50', 'Median Length']
         existing_asm = [c for c in cols_asm if c in clean_mdbg.columns]
         clean_mdbg = clean_mdbg[existing_asm]
         
         # EXPORT
         mdbg_csv = os.path.join(out_dir, 'assembly_stats_mdbg.csv')
         clean_mdbg.to_csv(mdbg_csv, index=False)
         print(f"Exported MetaMDBG stats to: {mdbg_csv}")

         html_body += f'<div class="table-wrapper">{df_to_html_table(clean_mdbg, "mdbg-table")}</div>'
    else:
         html_body += "<p>No MetaMDBG stats available.</p>"
else:
    html_body += "<p>No MetaMDBG stats available.</p>"
html_body += '</div>'

# Flye Table
html_body += '<div class="tab-content" id="stats-flye">'
if df_stats_flye is not None:
    clean_flye = clean_assembly_stats(df_stats_flye) if df_stats_flye is not None else None
    if clean_flye is not None:
         # Keep Contigs, N50, Median Length
         cols_asm = ['Sample', 'Contigs', 'N50', 'Median Length']
         existing_asm = [c for c in cols_asm if c in clean_flye.columns]
         clean_flye = clean_flye[existing_asm]
         
         # EXPORT
         flye_csv = os.path.join(out_dir, 'assembly_stats_flye.csv')
         clean_flye.to_csv(flye_csv, index=False)
         print(f"Exported Flye stats to: {flye_csv}")

         html_body += f'<div class="table-wrapper">{df_to_html_table(clean_flye, "flye-table")}</div>'
    else:
         html_body += "<p>No Flye stats available.</p>"
else:
    html_body += "<p>No Flye stats available.</p>"
html_body += '</div>'

# Myloasm Table
html_body += '<div class="tab-content" id="stats-myloasm">'
if df_stats_myloasm is not None:
    clean_myloasm = clean_assembly_stats(df_stats_myloasm) if df_stats_myloasm is not None else None
    if clean_myloasm is not None:
         cols_asm = ['Sample', 'Contigs', 'N50', 'Median Length']
         existing_asm = [c for c in cols_asm if c in clean_myloasm.columns]
         clean_myloasm = clean_myloasm[existing_asm]
         
         # EXPORT
         myloasm_csv = os.path.join(out_dir, 'assembly_stats_myloasm.csv')
         clean_myloasm.to_csv(myloasm_csv, index=False)
         print(f"Exported Myloasm stats to: {myloasm_csv}")

         html_body += f'<div class="table-wrapper">{df_to_html_table(clean_myloasm, "myloasm-table")}</div>'
    else:
         html_body += "<p>No Myloasm stats available.</p>"
else:
    html_body += "<p>No Myloasm stats available.</p>"
html_body += '</div></div>'


# 5. Listeria Contigs (New Section)
html_body += '<div class="section"><h2>Listeria Contigs</h2>'
if df_listeria is not None and len(df_listeria) > 0 and ('contigs_flye' in df_listeria.columns or 'contigs_mdbg' in df_listeria.columns or 'contigs_myloasm' in df_listeria.columns):
    html_body += '<h3>Listeria Contigs Data (Flye, MetaMDBG & Myloasm)</h3>'
    # Display Columns: Sample, Type, Listeria Contigs (Flye), Listeria Contigs (MDBG)
    # Plus new metric columns: bases_flye, median_flye, bases_mdbg, median_mdbg
    
    display_cols_contigs = ['sample', 'type', 
                            'Flye Contigs', 'flye_contig_bases', 'flye_median_contig_len',
                            'MetaMDBG Contigs', 'mdbg_contig_bases', 'mdbg_median_contig_len',
                            'Myloasm Contigs', 'myloasm_contig_bases', 'myloasm_median_contig_len']
    existing_cols_contigs = [c for c in display_cols_contigs if c in df_listeria.columns]
    
    final_df_contigs = df_listeria[existing_cols_contigs].rename(columns={
        'sample': 'Sample',
        'type': 'Type',
        'MetaMDBG Contigs': 'MDBG Contigs',
        'flye_contig_bases': 'Flye Bases',
        'flye_median_contig_len': 'Flye Median',
        'mdbg_contig_bases': 'MDBG Bases',
        'mdbg_median_contig_len': 'MDBG Median',
        'myloasm_contig_bases': 'Myloasm Bases',
        'myloasm_median_contig_len': 'Myloasm Median'
    })
    
    # FILTER: Remove rows where ALL are 0
    cond_cols = [c for c in ['Flye Contigs', 'MDBG Contigs', 'Myloasm Contigs'] if c in final_df_contigs.columns]
    if cond_cols:
        final_df_contigs = final_df_contigs[final_df_contigs[cond_cols].sum(axis=1) > 0].copy()
    
    # EXPORT
    listeria_contigs_csv = os.path.join(out_dir, 'listeria_contigs_filtered.csv')
    final_df_contigs.to_csv(listeria_contigs_csv, index=False)
    print(f"Exported Listeria Contigs table to: {listeria_contigs_csv}")

    html_body += f'<div class="search-box"><input type="text" placeholder="Search samples..." onkeyup="filterTable(this, \'listeria-contigs-table\')"></div>'
    html_body += f'<div class="table-wrapper">{df_to_html_table(final_df_contigs, "listeria-contigs-table")}</div>'
else:
    html_body += "<p>No Listeria contig counts available (run Step 14).</p>"
html_body += '</div>'


# 6. AMR Section (Last)
html_body += """
<div class="section">
  <h2>AMR Gene Detection</h2>
  
  <div class="tabs">
    <button class="tab-btn active" onclick="openTab(event, 'amr-reads')">Genes in Reads</button>
    <button class="tab-btn" onclick="openTab(event, 'amr-contigs')">Genes in Contigs</button>
  </div>
"""

# AMR from Reads
html_body += '<div class="tab-content active" id="amr-reads">'
if df_amr_reads is not None:
    # --- AMR Comparisons Plot (Reads) ---
    # Count unique genes per type
    if 'type' not in df_amr_reads.columns and 'Sample' in df_amr_reads.columns:
        df_amr_reads['type'] = df_amr_reads['Sample'].apply(lambda x: 'AS' if '_AS' in x else 'N')
            
    if 'Gene Symbol' in df_amr_reads.columns and 'type' in df_amr_reads.columns:
        # Group by type, nunique gene symbol
        amr_counts_reads = df_amr_reads.groupby('type')['Gene Symbol'].nunique().reset_index()
        amr_counts_reads.columns = ['type', 'Gene Count']
        
        p_amr_reads = create_vertical_bar_chart(amr_counts_reads, 'type', 'Gene Count', 'AMR Genes (Reads) AS vs N', 'Count')
        html_body += f'<div style="max-width:600px; margin-bottom:2rem;"><img src="data:image/png;base64,{p_amr_reads}" style="width:100%; border:1px solid #e2e8f0; border-radius:4px;"></div>'
        
    # EXPORT
    amr_reads_out = os.path.join(out_dir, 'amr_genes_reads.csv')
    df_amr_reads.to_csv(amr_reads_out, index=False)
    print(f"Exported AMR Reads table to: {amr_reads_out}")

    html_body += f'<div class="table-wrapper">{df_to_html_table(df_amr_reads, "amr-reads-table", max_rows=500)}</div>'
else:
    html_body += "<p>No AMR genes found in reads.</p>"
html_body += '</div>'

# AMR from Contigs
html_body += '<div class="tab-content" id="amr-contigs">'
if df_amr_contigs is not None:
    # --- AMR Comparisons Plot (Contigs) ---
    if 'type' not in df_amr_contigs.columns and 'Sample' in df_amr_contigs.columns:
        df_amr_contigs['type'] = df_amr_contigs['Sample'].apply(lambda x: 'AS' if '_AS' in x else 'N')
        
    if 'Gene Symbol' in df_amr_contigs.columns and 'type' in df_amr_contigs.columns:
         amr_counts_contigs = df_amr_contigs.groupby('type')['Gene Symbol'].nunique().reset_index()
         amr_counts_contigs.columns = ['type', 'Gene Count']
         p_amr_contigs = create_vertical_bar_chart(amr_counts_contigs, 'type', 'Gene Count', 'AMR Genes (Contigs) AS vs N', 'Count')
         html_body += f'<div style="max-width:600px; margin-bottom:2rem;"><img src="data:image/png;base64,{p_amr_contigs}" style="width:100%; border:1px solid #e2e8f0; border-radius:4px;"></div>'

    # EXPORT
    amr_contigs_out = os.path.join(out_dir, 'amr_genes_contigs.csv')
    df_amr_contigs.to_csv(amr_contigs_out, index=False)
    print(f"Exported AMR Contigs table to: {amr_contigs_out}")

    html_body += f'<div class="table-wrapper">{df_to_html_table(df_amr_contigs, "amr-contigs-table")}</div>'
else:
    html_body += "<p>No AMR genes found in contigs.</p>"
html_body += '</div></div>'


# Footer & Scripts
html_footer = f"""
<div class="footer">
  Raw data processed by Tim Thilo Maria Reska's Pipeline v1.0.<br>
  Report generated via customized Python engine.
</div>

</div> <!-- End Container -->

<script>
function openTab(evt, tabName) {{
  var i, tabcontent, tablinks;
  tabcontent = document.getElementsByClassName("tab-content");
  for (i = 0; i < tabcontent.length; i++) {{
    tabcontent[i].style.display = "none";
    tabcontent[i].classList.remove("active");
  }}
  tablinks = document.getElementsByClassName("tab-btn");
  for (i = 0; i < tablinks.length; i++) {{
    tablinks[i].className = tablinks[i].className.replace(" active", "");
  }}
  document.getElementById(tabName).style.display = "block";
  document.getElementById(tabName).classList.add("active");
  evt.currentTarget.className += " active";
}}

function sortTable(th) {{
  var table = th.closest("table");
  var tbody = table.querySelector("tbody");
  var rows = Array.from(tbody.querySelectorAll("tr"));
  var index = Array.from(th.parentNode.children).indexOf(th);
  var asc = th.getAttribute("data-order") !== "asc";
  
  rows.sort(function(a, b) {{
    var valA = a.children[index].innerText;
    var valB = b.children[index].innerText;
    
    var numA = parseFloat(valA.replace(/,/g, ""));
    var numB = parseFloat(valB.replace(/,/g, ""));
    
    if (!isNaN(numA) && !isNaN(numB)) {{
      return asc ? numA - numB : numB - numA;
    }}
    return asc ? valA.localeCompare(valB) : valB.localeCompare(valA);
  }});
  
  rows.forEach(function(row) {{ tbody.appendChild(row); }});
  th.setAttribute("data-order", asc ? "asc" : "desc");
}}

function filterTable(input, tableId) {{
  var filter = input.value.toUpperCase();
  var table = document.getElementById(tableId);
  var tr = table.getElementsByTagName("tr");
  for (var i = 1; i < tr.length; i++) {{
    var td = tr[i].getElementsByTagName("td")[0]; // Filter by first column (Sample)
    if (td) {{
      var txtValue = td.textContent || td.innerText;
      if (txtValue.toUpperCase().indexOf(filter) > -1) {{
        tr[i].style.display = "";
      }} else {{
        tr[i].style.display = "none";
      }}
    }}
  }}
}}
</script>
</body>
</html>
"""

# Assemble Full Report
full_report = html_head + html_body + html_footer

# Assemble Slim Report (Remove Methods Section)
# Methods section is wrapped in <!-- BEGIN_METHODS --> ... <!-- END_METHODS -->
# Split by marker
slim_body = html_body.split('<!-- BEGIN_METHODS -->')[0] + html_body.split('<!-- END_METHODS -->')[1]
slim_report = html_head + slim_body + html_footer

# Write Output
out_file = os.path.join(out_dir, 'pipeline_report.html')
with open(out_file, 'w') as f:
    f.write(full_report)

out_file_slim = os.path.join(out_dir, 'pipeline_report_slim.html')
with open(out_file_slim, 'w') as f:
    f.write(slim_report)

print(f"Report generated: {out_file}")
print(f"Slim report generated: {out_file_slim}")

# Generate PDF (Slim only? or Full? User asked for 'all the information', so Full)
generate_pdf(out_file)
generate_pdf(out_file_slim)

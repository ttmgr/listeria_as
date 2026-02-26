#!/usr/bin/env python3
"""
Step 19 reads-focused report.
Purpose: create a quick HTML summary of read-level metrics before assembly completes.

Outputs:
  processing/report/reads_report.html
  processing/report/reads_methods.md

Usage:
  python3 19_reads_report.py <base_dir>
"""

import sys
import os
import glob
import subprocess
import pandas as pd
import numpy as np
from datetime import datetime

base_dir = sys.argv[1]
out_dir = os.path.join(base_dir, 'processing/report')
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
                for line in out.split('\n'):
                    line = line.strip()
                    if line and any(c.isdigit() for c in line):
                        return line
                return out.split('\n')[0]
        except Exception:
            continue
    return 'N/A'

def df_to_html_table(df, table_id='', max_rows=None):
    if df is None or df.empty:
        return '<p class="text-muted">No data available</p>'
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
# Collect tool versions (Reads steps only)
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
    ('Kraken2', 'kraken2', 'Taxonomic classification (reads)', '05',
     'kraken2 --db DB --use-names --threads 20 --report report.txt --output output.txt input.fastq'),
    ('seqtk', 'seqtk', 'Listeria read extraction', '06',
     'seqtk subseq input.fastq read_ids.txt > listeria.fastq'),
    ('seqkit', 'seqkit', 'FASTQ → FASTA conversion', '10',
     'seqkit fq2fa input.fastq > output.fasta'),
    ('AMRFinderPlus', 'amrfinder', 'AMR gene detection (reads)', '11',
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

# 1. Read Metrics
read_csv = os.path.join(base_dir, 'processing/stats/read_metrics_summary.csv')
df_reads = pd.read_csv(read_csv) if os.path.exists(read_csv) else pd.DataFrame()

# 2. Listeria Reads (Reads Only)
# Try main overview first, else fallback to summary
listeria_csv = os.path.join(base_dir, 'processing/listeria/overview/listeria_overview.csv')
df_listeria = None
if os.path.exists(listeria_csv):
    df_temp = pd.read_csv(listeria_csv)
    # Filter only relevant read columns
    cols = ['Sample', 'Listeria Reads', 'Listeria (%)', 'Total Reads', 'Barcode', 'Type']
    # Check if we have polished headers or old ones
    if 'Listeria Reads' in df_temp.columns:
         df_listeria = df_temp[cols].copy()
    else:
        # Map old to new
        df_listeria = df_temp[['sample', 'listeria_reads', 'pct_listeria', 'total_reads', 'barcode', 'type']].rename(columns={
            'sample': 'Sample',
            'listeria_reads': 'Listeria Reads',
            'pct_listeria': 'Listeria (%)',
            'total_reads': 'Total Reads',
            'barcode': 'Barcode',
            'type': 'Type'
        })
else:
    # Fallback to listeria_summary.tsv (raw)
    listeria_tsv = os.path.join(base_dir, 'processing/listeria/listeria_summary.tsv')
    if os.path.exists(listeria_tsv):
        df_listeria = pd.read_csv(listeria_tsv, sep='\t', header=None,
                                  names=['Sample', 'Listeria Reads', 'Listeria Bases', 'Mean Len', 'Median Len'])
        # Drop duplicates
        df_listeria = df_listeria.drop_duplicates(subset='Sample', keep='last')

# 3. AMR Reads
# Try overview csv
amr_reads_csv = os.path.join(base_dir, 'processing/amrfinder/overview/amr_reads_overview.csv')
df_amr = None
if os.path.exists(amr_reads_csv):
    df_amr = pd.read_csv(amr_reads_csv)
else:
    # Try parsing raw files if summary missing
    pass # Too complex to duplicate logic here, assume step 16 runs or user waits

# ============================================================
# Build HTML
# ============================================================

print("Generating HTML report...")

# KPIs
n_samples = len(df_reads)
total_reads = df_reads['number_of_reads'].sum() if not df_reads.empty else 0
listeria_hits = df_listeria['Listeria Reads'].sum() if df_listeria is not None else 0
amr_hits = df_amr['Read Count'].sum() if df_amr is not None and 'Read Count' in df_amr.columns else 0

css = """
  :root { --bg: #ffffff; --text: #18181b; --accent: #2563eb; --border: #e4e4e7; --muted: #71717a; }
  * { box-sizing: border-box; }
  body { font-family: -apple-system, sans-serif; background: #f8fafc; color: var(--text); padding: 2rem; line-height: 1.6; }
  .container { max-width: 1200px; margin: 0 auto; background: white; padding: 2rem; border-radius: 8px; box-shadow: 0 1px 3px rgba(0,0,0,0.1); }
  h1 { font-size: 1.8rem; margin-bottom: 0.5rem; color: #0f172a; }
  h2 { font-size: 1.4rem; margin-top: 2rem; margin-bottom: 1rem; border-bottom: 2px solid var(--border); padding-bottom: 0.5rem; }
  .subtitle { color: var(--muted); margin-bottom: 2rem; }
  
  .kpi-grid { display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr)); gap: 1rem; margin-bottom: 2rem; }
  .kpi-card { background: #f1f5f9; padding: 1rem; border-radius: 6px; text-align: center; }
  .kpi-card .value { font-size: 1.8rem; font-weight: 700; color: var(--accent); }
  .kpi-card .label { font-size: 0.85rem; color: var(--muted); text-transform: uppercase; }

  table { width: 100%; border-collapse: collapse; font-size: 0.85rem; margin-top: 0.5rem; }
  th { background: #f8fafc; text-align: left; padding: 0.6rem; border-bottom: 2px solid var(--border); cursor: pointer; }
  td { padding: 0.6rem; border-bottom: 1px solid var(--border); }
  tr:hover { background: #f8fafc; }
  .table-wrapper { overflow-x: auto; border: 1px solid var(--border); border-radius: 6px; max-height: 500px; }
  
  .footer { margin-top: 3rem; pt: 1rem; border-top: 1px solid var(--border); font-size: 0.8rem; color: var(--muted); text-align: center; }
"""

html = f"""<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">
<title>Nanopore Reads Report</title>
<style>{css}</style>
</head>
<body>
<div class="container">

<div class="header">
  <h1>Nanopore Reads QC Report</h1>
  <div class="subtitle">Generated: {datetime.now().strftime('%Y-%m-%d %H:%M')}</div>
</div>

<div class="kpi-grid">
  <div class="kpi-card"><div class="value">{n_samples}</div><div class="label">Samples</div></div>
  <div class="kpi-card"><div class="value">{total_reads:,.0f}</div><div class="label">Total Reads</div></div>
  <div class="kpi-card"><div class="value">{listeria_hits:,.0f}</div><div class="label">Listeria Reads</div></div>
  <div class="kpi-card"><div class="value">{amr_hits:,.0f}</div><div class="label">AMR Read Hits</div></div>
</div>

<div class="section">
  <h2>Processing Methods</h2>
  <div class="table-wrapper">
    {df_to_html_table(pd.DataFrame(tool_versions))}
  </div>
</div>

<div class="section">
  <h2>Read Statistics (QC)</h2>
  <div class="table-wrapper">
    {df_to_html_table(df_reads, 'metrics-table')}
  </div>
</div>

<div class="section">
  <h2>Listeria Detections (Reads)</h2>
  <div class="table-wrapper">
    {df_to_html_table(df_listeria, 'listeria-table')}
  </div>
</div>

<div class="section">
  <h2>AMR Gene Detections (Reads)</h2>
  <div class="table-wrapper">
    {df_to_html_table(df_amr, 'amr-table')}
  </div>
</div>

<div class="footer">
  Reads-Only Report | Urban Lab, Helmholtz Munich
</div>

</div>
<script>
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
</html>
"""

with open(os.path.join(out_dir, 'reads_report.html'), 'w') as f:
    f.write(html)
print(f"Saved: {os.path.join(out_dir, 'reads_report.html')}")

# Save methods markdown
with open(os.path.join(out_dir, 'reads_methods.md'), 'w') as f:
    f.write("# Nanopore Read Processing Methods\n\n")
    f.write(pd.DataFrame(tool_versions).to_markdown(index=False))
print(f"Saved: {os.path.join(out_dir, 'reads_methods.md')}")

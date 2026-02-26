#!/usr/bin/env python3
"""
Step 22: Local Visualization — PhD-Quality Figures & Summary Tables
Generates SEPARATE reports for Black and Blue colour groups.

Usage:
    python3 22_local_plots.py <base_dir>
"""

import sys, os, re, io, base64
from datetime import datetime
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.patches import Patch
from matplotlib.backends.backend_pdf import PdfPages

# ============================================================
# CONFIG
# ============================================================

COLOUR_GROUPS = {
    'Black': {
        'title': 'Black Samples',
        'barcodes': [3,4,5,6,7, 14,15,16,18,19, 26,27,28,29,30],
        'sample_map': {3:'A1',4:'A2',5:'A3',6:'A4',7:'A5',
                       14:'C1',15:'C2',16:'C3',18:'C4',19:'C5',
                       26:'D1',27:'D2',28:'D3',29:'D4',30:'D5'},
        'outdir': 'plots_black',
    },
    'Blue': {
        'title': 'Blue Samples',
        'barcodes': [1,2,9,17,25, 8,10,11,12,13, 20,21,22,23,24],
        'sample_map': {1:'A1',2:'A5',9:'A2',17:'A3',25:'A4',
                       8:'C1',10:'C2',11:'C3',12:'C4',13:'C5',
                       20:'D1',21:'D2',22:'D3',23:'D4',24:'D5'},
        'outdir': 'plots_blue',
    },
}
LM_BARCODES = [31,32,33]
LM_MAP = {31:'Lm2',32:'Lm4',33:'Lm6'}

METHOD_MAP = {
    1:'A',2:'A',3:'A',4:'A',5:'A',6:'A',7:'A',9:'A',17:'A',25:'A',
    8:'C',10:'C',11:'C',12:'C',13:'C',14:'C',15:'C',16:'C',18:'C',19:'C',
    20:'D',21:'D',22:'D',23:'D',24:'D',26:'D',27:'D',28:'D',29:'D',30:'D',
    31:'Lm',32:'Lm',33:'Lm',
}

METHODS = ['A','C','D']
METHOD_FULL = {'A':'Sponge / PowerSoil','C':'Cotton / Zymo','D':'Zymo swab / Zymo','Lm':'Lm control / Zymo'}

# N = solid (velvet red), AS = hatched (velvet green)
C_N       = '#b03060'
C_N_FILL  = '#b0306040'
C_AS      = '#2d8a4e'
C_AS_FILL = '#2d8a4e40'

GRP_C   = {'A':'#3366cc','C':'#cc6633','D':'#339966','Lm':'#9966cc'}
GRP_C_L = {'A':'#3366cc30','C':'#cc663330','D':'#33996630','Lm':'#9966cc30'}

plt.rcParams.update({
    'font.family':'sans-serif','font.sans-serif':['Arial','Helvetica','DejaVu Sans'],
    'font.size':9,'axes.titlesize':10,'axes.titleweight':'bold','axes.labelsize':9,
    'axes.linewidth':0.7,'axes.spines.top':False,'axes.spines.right':False,
    'xtick.labelsize':8,'ytick.labelsize':8,'xtick.major.width':0.5,'ytick.major.width':0.5,
    'xtick.major.size':3,'ytick.major.size':3,'legend.fontsize':8,'legend.frameon':False,
    'figure.dpi':300,'savefig.dpi':300,'savefig.bbox':'tight','savefig.pad_inches':0.08,
})

# ============================================================
# Helpers
# ============================================================

def bc_num(name):
    m = re.search(r'barcode(\d+)', str(name))
    return int(m.group(1)) if m else 0

def cond(name):
    if '_AS' in str(name): return 'AS'
    if '_N' in str(name): return 'N'
    return 'unknown'

def panel(ax, letter, x=-0.10, y=1.06):
    ax.text(x, y, letter, transform=ax.transAxes, fontsize=12, fontweight='bold', va='top', ha='left')

def fig_b64(fig):
    buf = io.BytesIO()
    fig.savefig(buf, format='png', bbox_inches='tight', dpi=200)
    plt.close(fig)
    buf.seek(0)
    return base64.b64encode(buf.read()).decode()

def try_load(bd, *paths, sep=',', hdr=True):
    for p in paths:
        fp = os.path.join(bd, p)
        if os.path.exists(fp):
            print(f"  ✓ {p}")
            kw = {'sep':sep, 'on_bad_lines':'skip'}
            if not hdr: kw['header'] = None
            return pd.read_csv(fp, **kw)
    print(f"  ✗ {paths[0]}")
    return None

def enrich_all(df, scol='Sample'):
    if df is None or len(df)==0: return None
    df = df.copy()
    df['bc'] = df[scol].apply(bc_num)
    df['cond'] = df[scol].apply(cond)
    df['method'] = df['bc'].map(METHOD_MAP)
    return df[(df['bc']>=1) & (df['bc']<=33)].copy()

def add_seps(ax, bcs, key_fn):
    prev = None
    for i, bc in enumerate(bcs):
        g = key_fn(bc)
        if g != prev and prev is not None:
            ax.axvline(x=i-0.5, color='#ccc', lw=0.5, zorder=0)
        prev = g

def draw_bars(ax, bcs, data_fn, w=0.35):
    x = np.arange(len(bcs))
    for i, bc in enumerate(bcs):
        grp = METHOD_MAP.get(bc, 'A')
        val_n  = data_fn(bc, 'N')
        val_as = data_fn(bc, 'AS')
        ax.bar(i - w/2, val_n, w, color=GRP_C[grp], edgecolor='white', linewidth=0.3)
        ax.bar(i + w/2, val_as, w, color=GRP_C_L[grp], edgecolor=GRP_C[grp], linewidth=0.3, hatch='///')
    return x

def std_legend():
    return [Patch(fc='#888', label='N (solid)'),
            Patch(fc='#ccc', ec='#888', hatch='///', label='AS (hatched)'),
            Patch(fc=GRP_C['A'], label=f'A — {METHOD_FULL["A"]}'),
            Patch(fc=GRP_C['C'], label=f'C — {METHOD_FULL["C"]}'),
            Patch(fc=GRP_C['D'], label=f'D — {METHOD_FULL["D"]}')]

# ============================================================
# Load data
# ============================================================

bd = sys.argv[1] if len(sys.argv) > 1 else '.'
print("Loading data...")

df_qc = try_load(bd, 'processing/report/qc_metrics.csv')
df_lr = try_load(bd, 'processing/report/listeria_reads_summary.csv')
df_af = try_load(bd, 'processing/report/assembly_stats_flye.csv')
df_am = try_load(bd, 'processing/report/assembly_stats_mdbg.csv')
df_ay = try_load(bd, 'processing/report/assembly_stats_myloasm.csv')

df_rl = try_load(bd, 'processing/read_lengths_filtered_agg.tsv', sep='\t', hdr=False)
if df_rl is not None:
    df_rl.columns = ['sample','length','state','count']
    df_rl = df_rl[df_rl['sample'].str.match(r'^barcode\d{2}_(AS|N)$', na=False)].copy()
    df_rl['length'] = pd.to_numeric(df_rl['length'], errors='coerce')
    df_rl['count']  = pd.to_numeric(df_rl['count'], errors='coerce')
    df_rl = df_rl.dropna(subset=['length','count'])
    df_rl['bc'] = df_rl['sample'].apply(bc_num)
    df_rl['cond'] = df_rl['sample'].apply(cond)
    df_rl['method'] = df_rl['bc'].map(METHOD_MAP)
    print(f"  → {len(df_rl)} clean rows")

df_k2 = try_load(bd, 'kraken2_reads.csv', 'processing/report/kraken2_reads.csv')

all_qc  = enrich_all(df_qc)
all_lr  = enrich_all(df_lr)

asm_all = []
for d, n in [(df_af,'Flye'),(df_am,'MetaMDBG'),(df_ay,'Myloasm')]:
    e = enrich_all(d)
    if e is not None:
        e['assembler'] = n
        asm_all.append(e)
all_asm = pd.concat(asm_all, ignore_index=True) if asm_all else None

# ============================================================
# GENERATE REPORTS PER COLOUR GROUP
# ============================================================

for colour, cfg in COLOUR_GROUPS.items():
    barcodes = cfg['barcodes']
    smap = cfg['sample_map']
    title = cfg['title']
    out = os.path.join(bd, cfg['outdir'])
    os.makedirs(out, exist_ok=True)

    # Include Lm controls alongside each group
    all_bcs = barcodes + LM_BARCODES
    full_smap = {**smap, **LM_MAP}

    def slabel(bc):
        return full_smap.get(bc, f'bc{bc:02d}')

    def sort_bcs(bcs_in):
        mo = {'A':0,'C':1,'D':2,'Lm':3}
        return sorted(bcs_in, key=lambda b: (mo.get(METHOD_MAP.get(b,'Lm'),9), b))

    # Filter data for this colour group + Lm
    qc  = all_qc[all_qc['bc'].isin(all_bcs)].copy() if all_qc is not None else None
    lr  = all_lr[all_lr['bc'].isin(all_bcs)].copy() if all_lr is not None else None
    asm = all_asm[all_asm['bc'].isin(all_bcs)].copy() if all_asm is not None else None
    rl  = df_rl[df_rl['bc'].isin(all_bcs)].copy() if df_rl is not None else None

    plots = {}
    saved = []

    print(f"\n{'='*60}")
    print(f"  GENERATING: {title} + Lm controls")
    print(f"  Output: {out}")
    print(f"{'='*60}")

    # ---- FIG 1: LISTERIA READS ----
    print(f"\n--- Fig 1: Listeria Reads ({colour}) ---")
    if lr is not None and len(lr) > 0:
        bcs = sort_bcs(lr['bc'].unique())
        fig, ax = plt.subplots(figsize=(8, 3.5))
        def lr_val(bc, cd):
            s = lr[(lr['bc']==bc) & (lr['cond']==cd)]
            return float(s['Listeria Reads'].values[0]) if len(s)>0 else 0
        x = draw_bars(ax, bcs, lr_val)
        add_seps(ax, bcs, lambda b: METHOD_MAP.get(b,''))
        ax.set_xticks(x)
        ax.set_xticklabels([slabel(b) for b in bcs], fontsize=7, rotation=45, ha='right')
        ax.set_ylabel('Listeria read count')
        ax.set_title(f'Listeria reads — {title}')
        ax.legend(handles=std_legend(), fontsize=6, loc='upper right', ncol=2)
        plt.tight_layout()
        p = os.path.join(out, 'fig1_listeria_reads.png')
        fig.savefig(p); fig.savefig(p.replace('.png','.pdf'))
        saved.append(p); plots['fig1'] = fig_b64(fig)
        print(f"  Saved: {os.path.basename(p)}")
    else:
        print("  SKIPPED")

    # ---- FIG 2: RELATIVE ABUNDANCE ----
    print(f"\n--- Fig 2: Relative Abundance ({colour}) ---")
    if lr is not None and len(lr) > 0:
        bcs = sort_bcs(lr['bc'].unique())
        fig, ax = plt.subplots(figsize=(8, 3.5))
        def pct_val(bc, cd):
            s = lr[(lr['bc']==bc) & (lr['cond']==cd)]
            return float(s['Percentage of Total Reads'].values[0]) if len(s)>0 else 0
        x = draw_bars(ax, bcs, pct_val)
        add_seps(ax, bcs, lambda b: METHOD_MAP.get(b,''))
        ax.set_xticks(x)
        ax.set_xticklabels([slabel(b) for b in bcs], fontsize=7, rotation=45, ha='right')
        ax.set_ylabel('Listeria reads (% of total)')
        ax.set_title(f'Relative abundance — {title}')
        ax.legend(handles=std_legend(), fontsize=6, loc='upper right', ncol=2)
        plt.tight_layout()
        p = os.path.join(out, 'fig2_relative_abundance.png')
        fig.savefig(p); fig.savefig(p.replace('.png','.pdf'))
        saved.append(p); plots['fig2'] = fig_b64(fig)
        print(f"  Saved: {os.path.basename(p)}")
    else:
        print("  SKIPPED")



    # ---- FIG 3: READ LENGTH DISTRIBUTIONS ----
    print(f"\n--- Fig 3: Read Lengths ({colour}) ---")
    if rl is not None and len(rl) > 0:
        maxl = int(rl['length'].max())
        bins = np.logspace(np.log10(50), np.log10(max(maxl, 100000)), 101)

        # Exclude Lm from main distributions
        rl_main = rl[rl['bc'].isin(barcodes)]

        # 3a: Combined AS vs N — red behind, green front
        fig, ax = plt.subplots(figsize=(7, 3.5))
        n_d  = rl_main[rl_main['cond']=='N']
        as_d = rl_main[rl_main['cond']=='AS']
        ax.hist(n_d['length'], bins=bins, weights=n_d['count'],
                histtype='stepfilled', color=C_N_FILL, edgecolor=C_N, linewidth=0.8,
                label='Normal', alpha=0.4, zorder=2)
        ax.hist(as_d['length'], bins=bins, weights=as_d['count'],
                histtype='stepfilled', color=C_AS_FILL, edgecolor=C_AS, linewidth=0.8,
                label='Adaptive Sampling', alpha=0.4, zorder=3)
        ax.axvline(x=400, color='#333', ls=':', lw=0.8, zorder=4)
        ax.set_xscale('log')
        ax.set_xlabel('Read length (bp)'); ax.set_ylabel('Read count')
        ax.set_title(f'Read length distribution — {title}')
        ax.set_xticks([100,1000,10000,100000]); ax.set_xticklabels(['100','1k','10k','100k'])
        ax.legend(fontsize=8)
        plt.tight_layout()
        p = os.path.join(out, 'fig3a_readlen_AS_vs_N.png')
        fig.savefig(p); fig.savefig(p.replace('.png','.pdf'))
        saved.append(p); plots['fig3a'] = fig_b64(fig)
        print(f"  Saved: {os.path.basename(p)}")

        # 3b: Per method — 3 panels
        fig, axes = plt.subplots(1, 3, figsize=(10, 3.2), sharey=True)
        for ax, grp in zip(axes, METHODS):
            g = rl_main[rl_main['method']==grp]
            n_g = g[g['cond']=='N']; as_g = g[g['cond']=='AS']
            ax.hist(n_g['length'], bins=bins, weights=n_g['count'],
                    histtype='stepfilled', color=C_N_FILL, edgecolor=C_N, linewidth=0.8,
                    label='N', alpha=0.4, zorder=2)
            ax.hist(as_g['length'], bins=bins, weights=as_g['count'],
                    histtype='stepfilled', color=C_AS_FILL, edgecolor=C_AS, linewidth=0.8,
                    label='AS', alpha=0.4, zorder=3)
            ax.axvline(x=400, color='#333', ls=':', lw=0.7)
            ax.set_xscale('log'); ax.set_xlabel('Read length (bp)')
            ax.set_title(f'{grp} — {METHOD_FULL[grp]}', fontsize=8, color=GRP_C[grp])
            ax.set_xticks([100,1000,10000,100000]); ax.set_xticklabels(['100','1k','10k','100k'])
            ax.legend(fontsize=6)
        axes[0].set_ylabel('Read count')
        for i,ax in enumerate(axes): panel(ax, chr(97+i))
        fig.suptitle(f'Read lengths per method — {title}', fontsize=10, fontweight='bold', y=1.04)
        plt.tight_layout()
        p = os.path.join(out, 'fig3b_readlen_per_method.png')
        fig.savefig(p); fig.savefig(p.replace('.png','.pdf'))
        saved.append(p); plots['fig3b'] = fig_b64(fig)
        print(f"  Saved: {os.path.basename(p)}")

        # 3c: Group comparison — transparent group colors
        fig, axes = plt.subplots(1, 2, figsize=(9, 3.5))
        for ax, cd, ttl in zip(axes, ['N','AS'], ['Normal','Adaptive Sampling']):
            for grp in METHODS:
                g = rl_main[(rl_main['method']==grp) & (rl_main['cond']==cd)]
                if len(g)==0: continue
                ax.hist(g['length'], bins=bins, weights=g['count'],
                        histtype='stepfilled', color=GRP_C_L[grp], edgecolor=GRP_C[grp],
                        linewidth=0.8, label=grp, alpha=0.4)
            ax.axvline(x=400, color='#333', ls=':', lw=0.7)
            ax.set_xscale('log'); ax.set_xlabel('Read length (bp)')
            ax.set_title(ttl, fontsize=9)
            ax.set_xticks([100,1000,10000,100000]); ax.set_xticklabels(['100','1k','10k','100k'])
            ax.legend(fontsize=7)
        axes[0].set_ylabel('Read count')
        panel(axes[0],'a'); panel(axes[1],'b')
        fig.suptitle(f'Method comparison — {title}', fontsize=10, fontweight='bold', y=1.04)
        plt.tight_layout()
        p = os.path.join(out, 'fig3c_readlen_groups.png')
        fig.savefig(p); fig.savefig(p.replace('.png','.pdf'))
        saved.append(p); plots['fig3c'] = fig_b64(fig)
        print(f"  Saved: {os.path.basename(p)}")

        # 3d: Lm controls separate
        rl_lm = rl[rl['bc'].isin(LM_BARCODES)]
        if len(rl_lm) > 0:
            fig, ax = plt.subplots(figsize=(7, 3.5))
            n_lm = rl_lm[rl_lm['cond']=='N']; as_lm = rl_lm[rl_lm['cond']=='AS']
            ax.hist(n_lm['length'], bins=bins, weights=n_lm['count'],
                    histtype='stepfilled', color=C_N_FILL, edgecolor=C_N, linewidth=0.8,
                    label='N', alpha=0.4, zorder=2)
            ax.hist(as_lm['length'], bins=bins, weights=as_lm['count'],
                    histtype='stepfilled', color=C_AS_FILL, edgecolor=C_AS, linewidth=0.8,
                    label='AS', alpha=0.4, zorder=3)
            ax.axvline(x=400, color='#333', ls=':', lw=0.7)
            ax.set_xscale('log'); ax.set_xlabel('Read length (bp)'); ax.set_ylabel('Read count')
            ax.set_title('Read lengths — Lm controls')
            ax.set_xticks([100,1000,10000,100000]); ax.set_xticklabels(['100','1k','10k','100k'])
            ax.legend(fontsize=8)
            plt.tight_layout()
            p = os.path.join(out, 'fig3d_readlen_Lm.png')
            fig.savefig(p); fig.savefig(p.replace('.png','.pdf'))
            saved.append(p); plots['fig3d'] = fig_b64(fig)
            print(f"  Saved: {os.path.basename(p)}")
    else:
        print("  SKIPPED")

    # ---- FIG 4: ASSEMBLY CONTIGS ----
    print(f"\n--- Fig 4: Assembly Contigs ({colour}) ---")
    if asm is not None and len(asm) > 0:
        assemblers = [a for a in ['Flye','MetaMDBG','Myloasm'] if a in asm['assembler'].values]
        na = len(assemblers)
        bcs = sort_bcs([b for b in asm['bc'].unique() if b in barcodes])
        # Add Lm separately at end
        lm_in_asm = [b for b in LM_BARCODES if b in asm['bc'].values]
        bcs_all = bcs + lm_in_asm
        fig, axes = plt.subplots(1, na, figsize=(3.5*na, 3.5), sharey=True)
        if na==1: axes = [axes]
        xp = np.arange(len(bcs_all))
        for idx,(ax,asmb) in enumerate(zip(axes, assemblers)):
            sub = asm[asm['assembler']==asmb]
            def asm_val(bc, cd, _s=sub):
                s = _s[(_s['bc']==bc)&(_s['cond']==cd)]
                return int(s['Contigs'].values[0]) if len(s)>0 else 0
            draw_bars(ax, bcs_all, asm_val)
            ax.set_xticks(xp)
            ax.set_xticklabels([slabel(b) for b in bcs_all], fontsize=6, rotation=45, ha='right')
            ax.set_title(asmb, fontweight='bold')
            panel(ax, chr(97+idx))
            add_seps(ax, bcs_all, lambda b: METHOD_MAP.get(b,''))
        axes[0].set_ylabel('Total contigs')
        fig.legend(handles=std_legend(), loc='lower center', ncol=5, fontsize=6, bbox_to_anchor=(0.5,-0.08))
        fig.suptitle(f'Assembly contigs — {title}', fontsize=10, fontweight='bold', y=1.02)
        plt.tight_layout()
        p = os.path.join(out, 'fig4_contigs.png')
        fig.savefig(p); fig.savefig(p.replace('.png','.pdf'))
        saved.append(p); plots['fig4'] = fig_b64(fig)
        print(f"  Saved: {os.path.basename(p)}")
    else:
        print("  SKIPPED")

    # ---- FIG 5: READ METRICS ----
    print(f"\n--- Fig 5: Read Metrics ({colour}) ---")
    if qc is not None and len(qc) > 0:
        mdefs = [(c,l) for c,l in [('Total Reads','Total reads'),('Median Length (bp)','Median length (bp)'),
                  ('Read N50 (bp)','Read N50 (bp)')] if c in qc.columns]
        nm = len(mdefs)
        bcs = sort_bcs([b for b in qc['bc'].unique() if b in barcodes])
        lm_qc = [b for b in LM_BARCODES if b in qc['bc'].values]
        bcs_all = bcs + lm_qc
        fig, axes = plt.subplots(1, nm, figsize=(3.5*nm, 3.5))
        if nm==1: axes = [axes]
        xp = np.arange(len(bcs_all))
        for idx,(ax,(col,ttl)) in enumerate(zip(axes, mdefs)):
            qc[col] = pd.to_numeric(qc[col], errors='coerce')
            def qc_val(bc, cd, _c=col):
                s = qc[(qc['bc']==bc)&(qc['cond']==cd)]
                return float(s[_c].values[0]) if len(s)>0 else 0
            draw_bars(ax, bcs_all, qc_val)
            ax.set_xticks(xp)
            ax.set_xticklabels([slabel(b) for b in bcs_all], fontsize=6, rotation=45, ha='right')
            ax.set_title(ttl, fontsize=9); panel(ax, chr(97+idx))
            if col=='Total Reads': ax.ticklabel_format(axis='y', style='sci', scilimits=(0,0))
            add_seps(ax, bcs_all, lambda b: METHOD_MAP.get(b,''))
        fig.suptitle(f'Read metrics — {title}', fontsize=10, fontweight='bold', y=1.02)
        fig.legend(handles=std_legend(), loc='lower center', ncol=5, fontsize=6, bbox_to_anchor=(0.5,-0.08))
        plt.tight_layout()
        p = os.path.join(out, 'fig5_read_metrics.png')
        fig.savefig(p); fig.savefig(p.replace('.png','.pdf'))
        saved.append(p); plots['fig5'] = fig_b64(fig)
        print(f"  Saved: {os.path.basename(p)}")
    else:
        print("  SKIPPED")

    # ---- KRAKEN2 TAXONOMY ----
    print(f"\n--- Kraken2 Taxonomy ({colour}) ---")
    if df_k2 is not None and len(df_k2) > 0:
        k2 = df_k2.copy()
        k2['bc'] = k2['barcode'].apply(bc_num)
        k2 = k2[k2['bc'].isin(all_bcs)].copy()
        k2['group'] = k2['bc'].map(METHOD_MAP)
        k2['method'] = k2['group'].map(METHOD_FULL)
        k2['sample'] = k2['bc'].apply(slabel)
        k2_out = k2.rename(columns={'type':'condition','taxon_name':'taxon','reads_direct':'reads'})
        k2_out = k2_out[['sample','group','method','condition','taxon','reads']].copy()
        k2_out['reads'] = pd.to_numeric(k2_out['reads'], errors='coerce').fillna(0).astype(int)
        k2_out = k2_out.sort_values(['group','sample','condition','reads'], ascending=[True,True,True,False])
        kp = os.path.join(out, 'kraken2_taxonomy.csv')
        k2_out.to_csv(kp, index=False)
        print(f"  Saved: {os.path.basename(kp)} ({len(k2_out)} rows)")
        k2_top = k2_out.groupby(['sample','condition']).apply(
            lambda x: x.nlargest(10, 'reads'), include_groups=False
        ).reset_index(drop=True)
        k2_top.to_csv(os.path.join(out, 'kraken2_top10.csv'), index=False)
    else:
        print("  SKIPPED")

    # ---- CSV SUMMARIES ----
    print(f"\n--- Summaries ({colour}) ---")
    rows = []
    for bc in sort_bcs(barcodes) + LM_BARCODES:
        for cd in ['N','AS']:
            row = {'Sample': slabel(bc), 'Group': METHOD_MAP.get(bc,''),
                   'Method': METHOD_FULL.get(METHOD_MAP.get(bc,''),''),
                   'Condition': cd, 'Barcode': bc}
            if qc is not None:
                q = qc[(qc['bc']==bc)&(qc['cond']==cd)]
                for c in ['Total Reads','Mean Length (bp)','Median Length (bp)','Read N50 (bp)','Total Bases (bp)']:
                    row[c] = q[c].values[0] if c in qc.columns and len(q)>0 else 0
            if lr is not None:
                l = lr[(lr['bc']==bc)&(lr['cond']==cd)]
                row['Listeria Reads'] = float(l['Listeria Reads'].values[0]) if len(l)>0 else 0
                row['% Listeria'] = float(l['Percentage of Total Reads'].values[0]) if len(l)>0 else 0
            if asm is not None:
                for asmb in ['Flye','MetaMDBG','Myloasm']:
                    a = asm[(asm['bc']==bc)&(asm['cond']==cd)&(asm['assembler']==asmb)]
                    row[f'{asmb} Contigs'] = int(a['Contigs'].values[0]) if len(a)>0 else 0
                    row[f'{asmb} N50'] = int(a['N50'].values[0]) if len(a)>0 else 0
            rows.append(row)
    df_sum = pd.DataFrame(rows)
    df_sum.to_csv(os.path.join(out, 'summary_per_sample.csv'), index=False)

    num_cols = [c for c in df_sum.columns if c not in ['Sample','Group','Method','Condition','Barcode']]
    # AS vs N (exclude Lm)
    agg_rows = []
    for cd in ['N','AS']:
        sub = df_sum[(df_sum['Condition']==cd)&(df_sum['Group']!='Lm')]
        row = {'Condition': cd, 'n': len(sub)}
        for c in num_cols:
            vals = pd.to_numeric(sub[c], errors='coerce').dropna()
            if len(vals)>0:
                row[f'{c} (median)'] = round(float(vals.median()),1)
                row[f'{c} (mean)'] = round(float(vals.mean()),1)
                row[f'{c} (SD)'] = round(float(vals.std()),1)
        agg_rows.append(row)
    df_agg = pd.DataFrame(agg_rows)
    df_agg.to_csv(os.path.join(out, 'summary_AS_vs_N.csv'), index=False)

    grp_rows = []
    for grp in METHODS + ['Lm']:
        for cd in ['N','AS']:
            sub = df_sum[(df_sum['Group']==grp)&(df_sum['Condition']==cd)]
            if len(sub)==0: continue
            row = {'Group': grp, 'Method': METHOD_FULL.get(grp,''), 'Condition': cd, 'n': len(sub)}
            for c in num_cols:
                vals = pd.to_numeric(sub[c], errors='coerce').dropna()
                if len(vals)>0:
                    row[f'{c} (median)'] = round(float(vals.median()),1)
                    row[f'{c} (mean)'] = round(float(vals.mean()),1)
            grp_rows.append(row)
    df_grp = pd.DataFrame(grp_rows)
    df_grp.to_csv(os.path.join(out, 'summary_per_group.csv'), index=False)
    print(f"  Saved 3 CSVs")

    # ---- COLOR-CODED PDF ----
    print(f"\n--- PDF Tables ({colour}) ---")
    ROW_COLORS = {'A':'#dbeafe','C':'#fef3c7','D':'#d1fae5','Lm':'#ede9fe'}
    try:
        pdf_path = os.path.join(out, 'summary_tables.pdf')
        with PdfPages(pdf_path) as pdf:
            dcols = [c for c in df_sum.columns if c != 'Barcode']
            fig_t, ax_t = plt.subplots(figsize=(16, 0.24*len(df_sum)+2))
            ax_t.axis('off')
            ax_t.set_title(f'Per-Sample Summary — {title}', fontsize=12, fontweight='bold', pad=12)
            ct, cc = [], []
            for _, row in df_sum[dcols].iterrows():
                ct.append([str(v) if not isinstance(v,float) else f'{v:,.1f}' for v in row])
                cc.append([ROW_COLORS.get(row.get('Group',''),'#fff')]*len(dcols))
            tbl = ax_t.table(cellText=ct, colLabels=dcols, cellColours=cc, loc='center', cellLoc='center')
            tbl.auto_set_font_size(False); tbl.set_fontsize(5)
            tbl.auto_set_column_width(list(range(len(dcols))))
            for c in range(len(dcols)):
                tbl[0,c].set_facecolor('#e0e0e0')
                tbl[0,c].set_text_props(fontweight='bold', fontsize=4.5)
            for (r,c),cell in tbl.get_celld().items():
                cell.set_edgecolor('#ccc'); cell.set_linewidth(0.3)
            plt.tight_layout(); pdf.savefig(fig_t); plt.close(fig_t)

            # AS vs N page
            cols2 = list(df_agg.columns)
            fig2, ax2 = plt.subplots(figsize=(14, 2.5)); ax2.axis('off')
            ax2.set_title(f'AS vs N — {title}', fontsize=12, fontweight='bold', pad=12)
            c2, cc2 = [], []
            for _, row in df_agg.iterrows():
                c2.append([str(v) if not isinstance(v,float) else f'{v:,.1f}' for v in row])
                cc2.append([('#e8f5e9' if row['Condition']=='AS' else '#fce4ec')]*len(cols2))
            t2 = ax2.table(cellText=c2, colLabels=cols2, cellColours=cc2, loc='center', cellLoc='center')
            t2.auto_set_font_size(False); t2.set_fontsize(5.5)
            t2.auto_set_column_width(list(range(len(cols2))))
            for c in range(len(cols2)):
                t2[0,c].set_facecolor('#e0e0e0'); t2[0,c].set_text_props(fontweight='bold', fontsize=5)
            for (r,c),cell in t2.get_celld().items():
                cell.set_edgecolor('#ccc'); cell.set_linewidth(0.3)
            plt.tight_layout(); pdf.savefig(fig2); plt.close(fig2)

        print(f"  Saved: {os.path.basename(pdf_path)}")
    except Exception as e:
        print(f"  PDF failed: {e}")

    # ---- HTML REPORT ----
    print(f"\n--- HTML Report ({colour}) ---")
    css = """
body{font-family:Arial,sans-serif;margin:0;background:#fafafa;color:#222}
.c{max-width:1100px;margin:0 auto;padding:24px;background:#fff}
.hd{background:#1a1a2e;color:white;padding:24px;border-radius:6px;margin-bottom:24px}
.hd h1{margin:0;font-size:1.3rem;letter-spacing:.5px}
.hd p{margin:6px 0 0;color:#8888aa;font-size:.85rem}
h2{color:#1a1a2e;border-bottom:1px solid #ddd;padding-bottom:6px;margin-top:28px;font-size:1.05rem}
img{max-width:100%;border:1px solid #eee;margin:12px 0}
.note{background:#f0f4ff;border-left:3px solid #4466cc;padding:10px 14px;font-size:.85rem;margin:12px 0}
pre{background:#f5f5f5;padding:12px;border-radius:4px;font-size:.8rem;overflow-x:auto}
.ft{text-align:center;padding:16px;color:#999;font-size:.75rem;border-top:1px solid #eee;margin-top:28px}
.lg{display:flex;gap:14px;flex-wrap:wrap;margin:8px 0;font-size:.8rem}
.lg span{display:flex;align-items:center;gap:4px}
.dot{width:12px;height:12px;border-radius:2px;display:inline-block}
"""
    html = f"""<!DOCTYPE html><html><head><meta charset="UTF-8">
<title>{title} — Extraction Method Comparison</title><style>{css}</style></head><body>
<div class="c">
<div class="hd">
<h1>{title} — Extraction Method Comparison</h1>
<p>Generated {datetime.now().strftime('%Y-%m-%d %H:%M')} · {colour} group + Lm controls</p>
</div>
<div class="lg">
<span><span class="dot" style="background:{GRP_C['A']}"></span> A — Sponge/PS</span>
<span><span class="dot" style="background:{GRP_C['C']}"></span> C — Cotton/Zymo</span>
<span><span class="dot" style="background:{GRP_C['D']}"></span> D — Zymo/Zymo</span>
<span><span class="dot" style="background:{GRP_C['Lm']}"></span> Lm control</span>
<span>| Solid=N, Hatched=AS</span>
</div>
"""
    for key, ftitle in [('fig1','1. Listeria Read Counts'),('fig2','2. Relative Abundance'),
                         ('fig3a','3a. Read Lengths — AS vs N'),('fig3b','3b. Read Lengths per Method'),
                         ('fig3c','3c. Method Comparison'),('fig3d','3d. Lm Controls'),
                         ('fig4','4. Assembly Contigs'),('fig5','5. Read Metrics')]:
        html += f'<h2>{ftitle}</h2>'
        if key in plots:
            html += f'<img src="data:image/png;base64,{plots[key]}" alt="{ftitle}">'
        else:
            html += '<div class="note">Data not available.</div>'

    html += f'<div class="ft">{datetime.now().strftime("%Y-%m-%d")} · {title}</div></div></body></html>'

    rp = os.path.join(out, 'report.html')
    with open(rp, 'w') as f: f.write(html)
    print(f"\n  HTML: {rp}")
    print(f"  Plots: {len(plots)}, Files: {len(saved)}")

print(f"\n{'='*60}")
print("ALL DONE!")
print(f"{'='*60}")

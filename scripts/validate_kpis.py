
import sys
import os
import pandas as pd
import numpy as np

def validate_kpis(base_dir):
    print(f"Validating data in: {base_dir}")
    
    # 1. Read Metrics (Total Reads, Total Bases)
    read_csv = os.path.join(base_dir, 'processing/stats/read_metrics_summary.csv')
    total_reads = 0
    total_bases_gb = 0
    if os.path.exists(read_csv):
        df_reads = pd.read_csv(read_csv)
        print("\n--- Read Metrics ---")
        if 'number_of_reads' in df_reads.columns:
            total_reads = df_reads['number_of_reads'].sum()
            print(f"Total Reads: {total_reads:,}")
        if 'total_bases' in df_reads.columns:
            total_bases_gb = df_reads['total_bases'].sum() / 1e9
            print(f"Total Bases: {total_bases_gb:.2f} Gb")
    else:
        print(f"\nError: {read_csv} not found.")

    # 2. Listeria Overview (Positives Breakdown)
    # The main script tries CSV first, then TSV. We should do the same to be robust.
    listeria_csv = os.path.join(base_dir, 'processing/listeria/overview/listeria_overview.csv')
    listeria_tsv = os.path.join(base_dir, 'processing/listeria/listeria_summary.tsv')
    df_listeria = None
    
    if os.path.exists(listeria_csv):
        print(f"\nLoading {listeria_csv}...")
        try:
            df_listeria = pd.read_csv(listeria_csv)
            # Rename if needed (as in main script)
            # The script logic:
            # rename_map = {'Sample': 'sample', 'Listeria Reads': 'listeria_reads', 'Listeria (%)': 'listeria_ratio', 'Type': 'type'}
            # Also handle old/raw names: {'pct_listeria': 'listeria_ratio', ...}
            
            # Simple normalization
            df_listeria.columns = df_listeria.columns.str.lower().str.replace(' ', '_').str.replace('(', '').str.replace(')', '').str.replace('%', 'ratio')
            # Check for 'listeria_reads'
            # If "Listeria Reads" -> "listeria_reads" (match)
            
            # If type not present, derive it
            if 'type' not in df_listeria.columns and 'sample' in df_listeria.columns:
                  df_listeria['type'] = df_listeria['sample'].apply(lambda x: 'AS' if '_AS' in x else 'N')
            
        except Exception as e:
            print(f"Error reading listeria_overview.csv: {e}")

    elif os.path.exists(listeria_tsv):
        print(f"\nLoading {listeria_tsv} (fallback)...")
        try:
             df_listeria = pd.read_csv(listeria_tsv, sep='\t', header=None)
             # Assume raw columns: sample, reads, bases, mean, median
             cols = ['sample', 'listeria_reads', 'listeria_bases', 'listeria_mean_len']
             if len(df_listeria.columns) == 5: cols.append('listeria_median_len')
             df_listeria.columns = cols[:len(df_listeria.columns)]
             df_listeria['type'] = df_listeria['sample'].apply(lambda x: 'AS' if '_AS' in x else 'N')
        except Exception as e:
            print(f"Error reading TSV: {e}")

    if df_listeria is not None:
        if 'listeria_reads' in df_listeria.columns:
            # Filter > 0
            positives = df_listeria[df_listeria['listeria_reads'] > 0]
            n_pos = len(positives)
            n_as = len(positives[positives['type'] == 'AS'])
            n_n = len(positives[positives['type'] == 'N'])
            
            print("--- Listeria Positives Breakdown ---")
            print(f"Total Positive Samples: {n_pos}")
            print(f"  - AS: {n_as}")
            print(f"  - N:  {n_n}")
        else:
            print("Error: 'listeria_reads' column missing.")
    else:
        print("\nError: Could not load Listeria Overview data.")

    # 3. Listeria Contigs (Filter Logic)
    contigs_tsv = os.path.join(base_dir, 'processing/listeria/listeria_contigs_summary.tsv')
    if os.path.exists(contigs_tsv):
        try:
            print(f"\nLoading {contigs_tsv}...")
            # Columns: BASENAME, ASSEMBLER, TOTAL_CONTIGS...
            df_lc = pd.read_csv(contigs_tsv, sep='\t', header=None, 
                                names=['sample', 'assembler', 'contigs', 'bases', 'median', 'k_count'])
            # Max to handle duplicates (same logic as script)
            df_lc = df_lc.groupby(['sample', 'assembler'], as_index=False).max()
            
            # Pivot
            df_lc_p = df_lc.pivot(index='sample', columns='assembler', values='contigs').reset_index()
            # Fillna 0
            df_lc_p = df_lc_p.fillna(0)
            
            print(f"Pivot columns (assemblers found): {df_lc_p.columns.tolist()}")
            
            # Check for flye / mdbg columns
            # The script uses exact names from the file. If file has 'flye' and 'metaMDBG', those are cols.
            flye_col = next((c for c in df_lc_p.columns if 'flye' in c.lower()), None)
            mdbg_col = next((c for c in df_lc_p.columns if 'meta' in c.lower() or 'mdbg' in c.lower()), None)
            
            if flye_col and mdbg_col:
                myloasm_col = next((c for c in df_lc_p.columns if 'myloasm' in c.lower()), None)
                # Logic: (Flye > 0) | (MetaMDBG > 0) | (Myloasm > 0)
                filter_expr = (df_lc_p[flye_col] > 0) | (df_lc_p[mdbg_col] > 0)
                filter_desc = f"({flye_col} > 0) | ({mdbg_col} > 0)"
                if myloasm_col:
                    filter_expr = filter_expr | (df_lc_p[myloasm_col] > 0)
                    filter_desc += f" | ({myloasm_col} > 0)"
                filtered = df_lc_p[filter_expr]
                print("\n--- Listeria Contigs Filter Logic ---")
                print(f"Total Samples (Pre-filter): {len(df_lc_p)}")
                print(f"Samples (Post-filter): {len(filtered)}")
                print(f"Filter used: {filter_desc}")
            else:
                 print(f"Warning: Could not identify Flye/MetaMDBG columns. Found: {df_lc_p.columns.tolist()}")

        except Exception as e:
            print(f"Error analyzing contigs: {e}")
    else:
        print(f"\nError: {contigs_tsv} not found.")

if __name__ == "__main__":
    if len(sys.argv) < 2:
        # Default to current dir's parent if not provided, assuming script is in preprocessing/
        # But user usually provides base dir
        print("Usage: python3 validate_kpis.py <base_dir>")
        sys.exit(1)
    base_dir = sys.argv[1]
    validate_kpis(base_dir)

# Adapting the Pipeline

The pipeline is highly modular and flexible. This guide covers how to prepare your metadata and how to alter the target extraction step to search for *any* organism, not just *Listeria monocytogenes*.

---

## 1. Using Sample Metadata (`sample_metadata.csv`)

For downstream reporting (like Python visualization models), you must connect abstract `barcode01` names to your corresponding sample attributes (kit, sample condition, DNA yield, etc.).

A template (`sample_metadata_template.csv`) is provided in the repository root.

**Example Structure:**
```csv
SampleName,Barcode,DNA_Concentration,Group,Swab,Kit,Cond,Type,TargetSpike
A1,barcode03,15.2,Black_A,Sponge,PowerSoil,AS,Sample,None
A1_N,barcode03,12.1,Black_A,Sponge,PowerSoil,N,Sample,None
```

### Applying Metadata in Python Scripts
To integrate the metadata, ensure your downstream visualization scripts (like `22_local_plots.py`) parse the CSV to map the generic output IDs to the `SampleName` column for axis labeling and comparative grouping.

---

## 2. Adapting for Other Organisms

This repository uses **Listeria** as its target for extraction (hence `listeria_as`). However, it is trivial to repurpose this for **Salmonella**, **E. coli**, **Campylobacter**, or any specific taxon identifiable by Kraken2.

### Modifying the Extraction Scripts
To target a different organism, you must adjust the search term regex logic in primarily two files. 

**File 1: `scripts/06_listeria_extract.sh`**
Look for the read extraction block reading the kraken outputs:
```bash
# Current
grep "Listeria" "$KRAKEN_DIR/${ID}_AS_kraken2_output.txt" | awk '{print $2}' > "$OUT_DIR/${ID}_AS_listeria_read_ids.txt"

# Change to target Salmonella (Note: capitalization is required!):
grep "Salmonella" "$KRAKEN_DIR/${ID}_AS_kraken2_output.txt" | awk '{print $2}' > "$OUT_DIR/${ID}_AS_salmonella_read_ids.txt"
```

**File 2: `scripts/14_listeria_contigs.sh`**
Similarly, for filtering assembled contigs:
```bash
# Current
grep "Listeria" "$KRAKEN_DIR/${ID}_AS_${asmb}_kraken2_output.txt" | awk '{print $2}' > "$OUT_DIR/${ID}_AS_${asmb}_listeria_contig_ids.txt"

# Change to target Pseudomonas (Note: capitalization is required!):
grep "Pseudomonas" "$KRAKEN_DIR/${ID}_AS_${asmb}_kraken2_output.txt" | awk '{print $2}' > "$OUT_DIR/${ID}_AS_${asmb}_pseudomonas_contig_ids.txt"
```

### Naming Conventions Downstream
While grep does the actual technical filtering, do not forget to rename the `listeria_overview.csv` summary compilation scripts (e.g. `15_compile_listeria_overview.sh`) if you wish your downstream reports to reflect the new target name.

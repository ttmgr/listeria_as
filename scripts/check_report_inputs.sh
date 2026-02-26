#!/bin/bash
# Check if all required inputs for the report exist and have content
# Usage: ./check_report_inputs.sh /lustre/groups/hpc/urban_lab/projects/tim

BASE_DIR="${1:-.}"
echo "Checking report inputs in: $BASE_DIR"
echo "--------------------------------------------------------"

check_file() {
    local f="$1"
    local desc="$2"
    if [[ -f "$f" ]]; then
        local size=$(du -h "$f" | cut -f1)
        local lines=$(wc -l < "$f")
        echo "[OK] $desc found ($size, $lines lines)"
        echo "     Path: $f"
        echo "     Header: $(head -n 1 "$f")"
    else
        echo "[MISSING] $desc NOT FOUND at: $f"
    fi
    echo ""
}

# 1. Read Metrics
check_file "$BASE_DIR/processing/stats/read_metrics_summary.csv" "Read Metrics Summary"

# 2. Listeria Overview
check_file "$BASE_DIR/processing/listeria/overview/listeria_overview.csv" "Listeria Overview CSV"

# 3. AMR Overview (Reads)
check_file "$BASE_DIR/processing/amrfinder/overview/amr_reads_overview.csv" "AMR Reads Overview"

# 4. AMR Overview (Contigs)
check_file "$BASE_DIR/processing/amrfinder/overview/amr_contigs_overview.csv" "AMR Contigs Overview"

# 5. Assembly Stats (MetaMDBG)
check_file "$BASE_DIR/processing/stats/assembly_stats_mdbg.tsv" "MetaMDBG Assembly Stats"

# 6. Assembly Stats (Flye)
check_file "$BASE_DIR/processing/stats/assembly_stats_flye.tsv" "Flye Assembly Stats"

# 7. Assembly Stats (Myloasm)
check_file "$BASE_DIR/processing/stats/assembly_stats_myloasm.tsv" "Myloasm Assembly Stats"

# 7. Listeria plots (check directory)
echo "Checking Listeria Plots directory:"
ls -lh "$BASE_DIR/processing/listeria/overview/" | grep ".png" | head -n 5
echo "(Total png files: $(ls -1 "$BASE_DIR/processing/listeria/overview/"*.png 2>/dev/null | wc -l))"
echo ""

echo "Done."

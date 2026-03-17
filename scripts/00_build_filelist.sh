#!/bin/bash
# -----------------------------------------------------------------------------
# Step 0: Build filelist.txt from round-based BAM directories.
#
# Scans listeria_1/, listeria_2/, … (as configured by ROUNDS in pipeline.conf)
# and writes one BAM path per line, relative to BAM_BASE_DIR.
#
# Usage: bash scripts/00_build_filelist.sh
# -----------------------------------------------------------------------------
SCRIPT_DIR="${SCRIPT_DIR:-$(cd "$(dirname "$0")" && pwd)}"
source "${SCRIPT_DIR}/pipeline.conf"

FILELIST_PATH="${WORK_DIR}/filelist.txt"
> "$FILELIST_PATH"   # truncate

echo "Building filelist from rounds: ${ROUNDS}"

for ROUND in $ROUNDS; do
    ROUND_DIR="${BAM_BASE_DIR}/listeria_${ROUND}"
    if [ ! -d "$ROUND_DIR" ]; then
        echo "WARNING: Round directory not found: $ROUND_DIR — skipping"
        continue
    fi
    # List all BAM files (exclude sorted BAMs), write paths relative to BAM_BASE_DIR
    find "$ROUND_DIR" -maxdepth 1 -type f -name '*.bam' ! -name '*.sorted.bam' | \
        sed "s#^${BAM_BASE_DIR}/##" | sort >> "$FILELIST_PATH"
done

NUM_FILES=$(wc -l < "$FILELIST_PATH")
echo "Filelist written to: $FILELIST_PATH"
echo "Total BAM files: $NUM_FILES"
echo ""
echo "First 5 entries:"
head -n 5 "$FILELIST_PATH"

# Verify BASENAME derivation for the first entry
if [ "$NUM_FILES" -gt 0 ]; then
    FIRST=$(head -n 1 "$FILELIST_PATH")
    echo ""
    echo "BASENAME derivation check:"
    echo "  Path:     $FIRST"
    echo "  BASENAME: $(derive_basename "$FIRST")"
fi

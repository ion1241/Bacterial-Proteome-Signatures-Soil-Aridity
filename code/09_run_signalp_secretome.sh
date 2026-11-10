#!/usr/bin/env bash
set -euo pipefail

mkdir -p signalp_out secretome transmem citoplasm

# Path to the SignalP 6.0 models
MODEL_DIR="/YOURDIRECTORY/signalp-6-package/models"

for faa in ./MAG_proteins/*.faa; do
    base="$(basename "$faa" .faa)"
    outdir="signalp_out/$base"

    mkdir -p "$outdir"

    echo "[INFO] Processing $faa ..."

    # Run SignalP 6.0
    OMP_NUM_THREADS=8 MKL_NUM_THREADS=8 \
    signalp6 \
        --fastafile "$faa" \
        --organism other \
        --mode fast \
        --format txt \
        --write_procs 8 \
        --bsize 64 \
        --output_dir "$outdir" \
        >"$outdir/run.log" 2>&1

    # Locate the prediction summary file
    sumfile="$outdir/prediction_results.txt"

    if [[ ! -s "$sumfile" ]]; then
        echo "[ERR] Failed for $base, see $outdir/run.log"
        continue
    fi

    # 1) Strict secretome (SP/Sec-SPI)
    awk -F'\t' '!/^#/ && $2=="SP"{print $1}' \
        "$sumfile" >> "$outdir/ids_secretome_strict.txt"

    # 2) Broadly secreted but non-SP proteins (LIPO, TAT, TATLIPO, PILIN)
    awk -F'\t' '!/^#/ && $2!="OTHER" && $2!="SP"{print $1}' \
        "$sumfile" >> "$outdir/ids_transmembrane.txt"

    # 3) Cytoplasmic proteins (OTHER)
    awk -F'\t' '!/^#/ && $2=="OTHER"{print $1}' \
        "$sumfile" >> "$outdir/ids_citoplasm.txt"

    # Split the FASTA according to the ID lists

    # 1) Strict secretome (SP/Sec-SPI)
    awk -v IDS="$outdir/ids_secretome_strict.txt" '
    BEGIN {
        while ((getline L < IDS) > 0) {
            sub(/\r$/, "", L)
            ids[L]=1
        }
    }
    /^>/ {
        hdr=$0
        sub(/^>/,"",hdr)
        sub(/[[:space:]]+$/,"",hdr)
        keep=(hdr in ids)
    }
    { if (keep) print }
    ' "$faa" >> "secretome/${base}.faa"

    # 2) Non-SP secreted proteins (LIPO/TAT/TATLIPO/PILIN)
    awk -v IDS="$outdir/ids_transmembrane.txt" '
    BEGIN {
        while ((getline L < IDS) > 0) {
            sub(/\r$/, "", L)
            ids[L]=1
        }
    }
    /^>/ {
        hdr=$0
        sub(/^>/,"",hdr)
        sub(/[[:space:]]+$/,"",hdr)
        keep=(hdr in ids)
    }
    { if (keep) print }
    ' "$faa" >> "transmem/${base}.faa"

    # 3) Cytoplasmic proteins (OTHER)
    awk -v IDS="$outdir/ids_citoplasm.txt" '
    BEGIN {
        while ((getline L < IDS) > 0) {
            sub(/\r$/, "", L)
            ids[L]=1
        }
    }
    /^>/ {
        hdr=$0
        sub(/^>/,"",hdr)
        sub(/[[:space:]]+$/,"",hdr)
        keep=(hdr in ids)
    }
    { if (keep) print }
    ' "$faa" >> "citoplasm/${base}.faa"

    echo "[OK] Generated:"
    echo " - secretome/${base}.faa"
    echo " - transmem/${base}.faa"
    echo " - citoplasm/${base}.faa"

done

echo "[DONE] Test completed."
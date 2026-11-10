#!/usr/bin/env python3
# -*- coding: utf-8 -*-

# Extracts and DECOMPRESSES .faa.gz files from gtdb_proteins_aa_reps_r226.tar.gz
# according to the accessions listed in ACCESSION_FILE, saving uncompressed
# .faa files in LOCAL_OUTPUT

import os
import re
import tarfile
import gzip
from pathlib import Path

# === PATHS: edit before running ===
TAR_PATH = Path("/YOURDIRECTORY/gtdb_proteins_aa_reps_r226.tar.gz")
LOCAL_OUTPUT = Path("/YOURDIRECTORY/faa")
ACCESSION_FILE = Path("/YOURDIRECTORY/soil_ids.txt")  # from 01_download_GTDB_soil_genomes.sh

# Create output directory if it does not already exist
os.makedirs(LOCAL_OUTPUT, exist_ok=True)

# Load accessions and normalize to GC[AF]_######### roots
# (remove RS_/GB_ prefixes and version suffixes such as .1/.2)
acc_re = re.compile(r'(?:RS_|GB_)?(GC[AF]_\d{9})(?:\.\d+)?')

roots = set()
with ACCESSION_FILE.open('r', encoding='utf-8') as f:
    for line in f:
        s = line.strip()
        if not s:
            continue
        m = acc_re.search(s)
        if m:
            roots.add(m.group(1))  # e.g. GCF_000195955

if not roots:
    raise SystemExit("No valid accessions found in selected_accessions.txt")

print(f"[INFO] Accessions loaded: {len(roots)}")

# Scan the tar archive and select entries containing one of the accession roots
# and ending in .faa.gz
files_to_extract = []

with tarfile.open(TAR_PATH, "r:gz") as tar:
    for member in tar.getmembers():
        if not member.isfile():
            continue
        if not member.name.endswith(".faa.gz"):
            continue

        # Example filename:
        # protein_faa_reps/bacteria/RS_GCF_000195955.2_protein.faa.gz
        m = acc_re.search(member.name)
        if m and m.group(1) in roots:
            files_to_extract.append(member)

print("Preview of files to extract:")
for m in files_to_extract[:10]:
    print(" ", m.name)

print(f"Total files to extract: {len(files_to_extract)}")

# === soil_ref_ids.txt: write normalized accessions without duplicates
mag_tsv_path = LOCAL_OUTPUT / "soil_ref_ids.txt"
written_accs = set()
mag_header_written = False

# Extract as uncompressed .faa files into LOCAL_OUTPUT
written = 0
skipped = 0

with tarfile.open(TAR_PATH, "r:gz") as tar, \
     mag_tsv_path.open("w", encoding="utf-8") as magfh:

    # Header
    magfh.write("accession\n")
    mag_header_written = True

    for member in files_to_extract:

        # Basename without directories and without .gz
        base_gz = Path(member.name).name   # RS_GCF_000195955.2_protein.faa.gz
        base_faa = re.sub(r'\.gz$', '', base_gz)  # RS_GCF_000195955.2_protein.faa
        out_path = LOCAL_OUTPUT / base_faa

        # Extract normalized accession from tar filename
        m_acc = acc_re.search(base_gz)
        if m_acc:
            acc_root = m_acc.group(1)  # GCF_######### / GCA_#########

            if acc_root not in written_accs:
                magfh.write(f"{acc_root}\n")
                written_accs.add(acc_root)

        # Avoid overwriting existing files
        if out_path.exists():
            skipped += 1
            continue

        # Read compressed bytes from the tar archive and decompress them
        fobj = tar.extractfile(member)
        if fobj is None:
            print(f"[WARN] Could not open {member.name} from the tar archive.")
            continue

        with gzip.GzipFile(fileobj=fobj, mode='rb') as gzf, \
             open(out_path, 'wb') as fout:

            # Copy in chunks
            while True:
                chunk = gzf.read(1024 * 1024)
                if not chunk:
                    break
                fout.write(chunk)

        written += 1

print(f"[INFO] Done. Wrote {written} decompressed .faa files to {LOCAL_OUTPUT}")

if skipped:
    print(f"[INFO] Skipped because they already existed: {skipped}")

# Quick listing
out_list = sorted(p.name for p in LOCAL_OUTPUT.glob("*.faa"))

print("Examples in output:")
for name in out_list[:10]:
    print(" ", name)

print("Total .faa files in output:", len(out_list))
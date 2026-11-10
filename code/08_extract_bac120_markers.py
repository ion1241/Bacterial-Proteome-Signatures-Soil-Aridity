import gzip, re
from pathlib import Path
from collections import defaultdict

# === PATHS: edit before running ===
# GTDB r226 bac120 marker genes (genomic_files_reps/bac120_marker_genes_reps_r226/faa/)
bac120_dir  = Path("/YOURDIRECTORY/bac120_marker_genes_reps_r226/faa/")
out_dir = Path("/YOURDIRECTORY/bac120_genomes/")
target_list = "/YOURDIRECTORY/target_mags.txt"  # one GTDB accession per line
out_dir.mkdir(parents=True, exist_ok=True)

def core_accession(acc):
    m = re.search(r'[GC][CA][AF]_(\d+)', acc)
    return m.group(1) if m else acc.strip()

with open(target_list) as f:
    raw_targets = [line.strip() for line in f if line.strip()]
target_cores = {core_accession(acc): acc for acc in raw_targets}
print(f"Target genomes: {len(target_cores)}", flush=True)

genome_seqs = defaultdict(list)
marker_files = sorted(bac120_dir.glob("*.faa.gz"))
print(f"Marker files found: {len(marker_files)}", flush=True)

for i, marker_faa in enumerate(marker_files, 1):
    marker_name = marker_faa.name.replace(".faa.gz", "")
    print(f"  [{i}/{len(marker_files)}] {marker_name}", flush=True)
    current_id, current_seq, current_core = None, [], None
    with gzip.open(marker_faa, "rt") as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                if current_core and current_seq and current_core in target_cores:
                    genome_seqs[current_id].append((marker_name, "".join(current_seq).rstrip("*")))
                current_id = line[1:].strip()
                current_core = core_accession(current_id)
                current_seq = []
            elif line:
                current_seq.append(line)
        if current_core and current_seq and current_core in target_cores:
            genome_seqs[current_id].append((marker_name, "".join(current_seq).rstrip("*")))

found_cores = {core_accession(gid) for gid in genome_seqs}
missing_cores = set(target_cores.keys()) - found_cores
print(f"\nMatched: {len(genome_seqs)} / {len(target_cores)} target genomes", flush=True)
if missing_cores:
    print(f"WARNING: {len(missing_cores)} targets not found in bac120 files:")
    for c in sorted(missing_cores):
        print(f"  {target_cores[c]} (core: {c})")

for genome_id, markers in genome_seqs.items():
    out_file = out_dir / f"{genome_id}_bac120.faa"
    with open(out_file, "w") as f:
        for marker_name, seq in markers:
            f.write(f">{genome_id}_{marker_name}\n{seq}\n")
    if len(markers) < 120:
        print(f"  WARNING: {genome_id} has only {len(markers)}/120 markers")

print(f"\nDone — {len(genome_seqs)} .faa files written to {out_dir}", flush=True)
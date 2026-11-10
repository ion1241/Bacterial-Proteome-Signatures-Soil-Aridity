#!/usr/bin/env python3

"""
per_protein_protparam (+ metapredict)
-------------------------------------
Processes a directory of FASTA files (.faa) and, for each file:

1. Computes per-contig metrics:
   - length
   - pI
   - acid_aa = percent_D + percent_E
   - basic_aa = percent_K + percent_H + percent_R

2. (Optional) Runs metapredict on the entire .faa file and computes
   six disorder-related variables per contig:
   - dis_mean = mean disorder score
   - dis_frac_0_5 = % of residues with score >= 0.5
   - dis_max_run_0_5 = longest continuous stretch with score >= 0.5
   - dis_top10_mean = mean of the top 10% most disordered residues
   - dis_N_term30 = mean disorder score over the first 30 residues
   - dis_C_term30 = mean disorder score over the last 30 residues

3. Writes one TSV per .faa file with rows = contigs and columns:
   seq_id, length, pI, acid_aa, basic_aa,
   dis_mean, dis_frac_0_5, dis_max_run_0_5,
   dis_top10_mean, dis_N_term30, dis_C_term30
"""

############################################################
# LIBRARIES
############################################################

import argparse                     # Read command-line arguments
import math
import multiprocessing as mp        # Parallel processing across faa files
from pathlib import Path            # Robust path handling
from typing import Dict, List

import pandas as pd                 # Build and save TSV files
from Bio import SeqIO               # Biopython FASTA reader
from Bio.SeqUtils.ProtParam import ProteinAnalysis
from tqdm import tqdm               # Progress bar
import sys                          # Write warnings to stderr
import metapredict as meta          # Per-residue disorder prediction

############################################################
# GENERAL CONFIGURATION
############################################################

STANDARD = list("ACDEFGHIKLMNPQRSTVWY")   # Standard amino acids
AMBIGUOUS = set("BJZXUO*-")               # Characters removed to avoid errors

############################################################
# 1. Sequence cleaning
############################################################

def clean_sequence(seq: str) -> str:
    """
    Receives a sequence string and:
      - Converts it to uppercase
      - Removes spaces and line breaks
      - Removes ambiguous characters

    Returns the cleaned sequence.
    """

    s = str(seq).upper().replace(" ", "").replace("\n", "")
    s = "".join(ch for ch in s if ch not in AMBIGUOUS)

    return s

############################################################
# 2. ProtParam per protein
############################################################

def per_protein_metrics(faa_path: Path) -> pd.DataFrame:
    """
    Reads a .faa file and returns a DataFrame with ONE ROW per protein.

    Columns:
      seq_id, length, pI, acid_aa, basic_aa
    """

    rows: List[Dict] = []

    # Iterate through all FASTA entries
    for rec in SeqIO.parse(str(faa_path), "fasta"):

        raw = str(rec.seq)
        aa = clean_sequence(raw)

        # Skip proteins that become empty after cleaning
        if not aa:
            continue

        row = {
            "seq_id": rec.id,
            "length": len(aa)
        }

        try:
            pa = ProteinAnalysis(aa)

            perc = pa.amino_acids_percent

            row.update({
                "pI": float(pa.isoelectric_point()),
                "acid_aa": float(
                    perc.get("D", 0.0) +
                    perc.get("E", 0.0)
                ),
                "basic_aa": float(
                    perc.get("K", 0.0) +
                    perc.get("H", 0.0) +
                    perc.get("R", 0.0)
                )
            })

        except Exception:
            row.update({
                "pI": math.nan,
                "acid_aa": math.nan,
                "basic_aa": math.nan
            })

        rows.append(row)

    return pd.DataFrame(rows)

############################################################
# 3. METAPREDICT
############################################################

def summarize_metapredict_scores(
    scores: List[float]
) -> Dict[str, float]:
    """
    Receives a list of disorder scores (one per residue, range 0–1)
    and calculates six disorder metrics.
    """

    L = len(scores)

    # If no scores are available, return an empty dictionary.
    # Missing values will later become NaN.
    if L == 0:
        return {}

    # Mean disorder level
    dis_mean = sum(scores) / L

    # Fraction of clearly disordered residues
    thr = 0.5
    dis_frac_0_5 = sum(
        1 for s in scores if s >= thr
    ) / L

    # Longest continuous disordered segment
    max_run = 0
    run = 0

    for s in scores:

        if s >= thr:
            run += 1

            if run > max_run:
                max_run = run

        else:
            run = 0

    # Mean disorder of the top 10% most disordered residues
    k = max(1, int(0.10 * L))

    top_scores = sorted(
        scores,
        reverse=True
    )[:k]

    dis_top10_mean = sum(top_scores) / len(top_scores)

    # N- and C-terminal disorder
    win = 30

    N_scores = scores[:win]
    C_scores = scores[-win:]

    dis_N_term30 = sum(N_scores) / len(N_scores)
    dis_C_term30 = sum(C_scores) / len(C_scores)

    return {
        "dis_mean": dis_mean,
        "dis_frac_0_5": dis_frac_0_5,
        "dis_max_run_0_5": max_run,
        "dis_top10_mean": dis_top10_mean,
        "dis_N_term30": dis_N_term30,
        "dis_C_term30": dis_C_term30,
    }


def run_metapredict_on_faa(
    faa_path: Path
) -> Dict[str, Dict[str, float]]:
    """
    Runs metapredict on the ENTIRE .faa file.

    Iterates through the output:

        header -> [sequence_str, scores_list]

    Returns:

        seq_id -> {dis_mean, dis_frac_0_5, ...}
    """

    pred_dict = meta.predict_disorder_fasta(
        str(faa_path)
    )

    out: Dict[str, Dict[str, float]] = {}

    for header, pair in pred_dict.items():

        # Keep only the sequence identifier from the FASTA header
        seq_id = header.split(" # ", 1)[0].strip()

        sequence_str, scores = pair

        out[seq_id] = summarize_metapredict_scores(
            scores
        )

    return out

############################################################
# 4. Worker: process ONE .faa file
############################################################

def process_file(args_tuple):

    (faa_path, out_dir, use_metapredict) = args_tuple

    # Convert strings to Path objects
    faa_path = Path(faa_path)

    out_tsv = out_dir / f"{faa_path.stem}.protparam.tsv"

    # 1) ProtParam → DataFrame
    df = per_protein_metrics(faa_path)

    # 2) metapredict (optional)
    if use_metapredict:

        try:

            met_map = run_metapredict_on_faa(faa_path)

            # Add six disorder columns by matching seq_id
            df["dis_mean"] = df["seq_id"].map(
                lambda s: met_map.get(
                    s, {}
                ).get("dis_mean", math.nan)
            )

            df["dis_frac_0_5"] = df["seq_id"].map(
                lambda s: met_map.get(
                    s, {}
                ).get("dis_frac_0_5", math.nan)
            )

            df["dis_max_run_0_5"] = df["seq_id"].map(
                lambda s: met_map.get(
                    s, {}
                ).get("dis_max_run_0_5", math.nan)
            )

            df["dis_top10_mean"] = df["seq_id"].map(
                lambda s: met_map.get(
                    s, {}
                ).get("dis_top10_mean", math.nan)
            )

            df["dis_N_term30"] = df["seq_id"].map(
                lambda s: met_map.get(
                    s, {}
                ).get("dis_N_term30", math.nan)
            )

            df["dis_C_term30"] = df["seq_id"].map(
                lambda s: met_map.get(
                    s, {}
                ).get("dis_C_term30", math.nan)
            )

        except Exception as e:

            # Continue even if metapredict fails
            print(
                f"[WARN] metapredict failed for {faa_path}: {e}",
                file=sys.stderr
            )

    # 3) Save final TSV
    df.to_csv(
        out_tsv,
        sep="\t",
        index=False
    )

    # Return summary information
    return {
        "file": str(faa_path),
        "out_tsv": str(out_tsv),
        "n_rows": len(df)
    }

############################################################
# 5. main(): prepare jobs and launch the pool
############################################################

def main():

    ap = argparse.ArgumentParser()

    ap.add_argument(
        "--input-dir",
        required=True,
        help="Directory containing .faa files"
    )

    ap.add_argument(
        "--glob",
        default="*.faa",
        help="File pattern (default: *.faa)"
    )

    ap.add_argument(
        "--out-dir",
        required=True,
        help="Output directory for TSV files"
    )

    ap.add_argument(
        "--workers",
        type=int,
        default=1,
        help="Number of parallel processes"
    )

    ap.add_argument(
        "--no-metapredict",
        action="store_true",
        help="Do not run metapredict (ProtParam only)"
    )

    args = ap.parse_args()

    in_dir = Path(args.input_dir)
    out_dir = Path(args.out_dir)

    out_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # List .faa files
    paths = sorted(
        in_dir.glob(args.glob)
    )

    if not paths:
        raise SystemExit(
            f"No files matched {in_dir}/{args.glob}"
        )

    use_metapredict = not args.no_metapredict

    # Prepare one job per .faa file
    jobs = [
        (p, out_dir, use_metapredict)
        for p in paths
    ]

    if args.workers > 1:

        # Parallel processing with multiprocessing.Pool
        with mp.Pool(
            processes=args.workers
        ) as pool:

            results = list(
                tqdm(
                    pool.imap(
                        process_file,
                        jobs
                    ),
                    total=len(jobs),
                    desc="FAA files"
                )
            )

    else:

        results = [
            process_file(j)
            for j in tqdm(
                jobs,
                desc="FAA files"
            )
        ]

    # Build a global index from worker outputs
    idx = pd.DataFrame(results)

    idx_path = out_dir / "_index.tsv"

    idx.to_csv(
        idx_path,
        sep="\t",
        index=False
    )

    print(
        f"[OK] Wrote {len(results)} TSVs to {out_dir}"
    )

    print(
        f"[OK] Index: {idx_path}"
    )


if __name__ == "__main__":
    main()
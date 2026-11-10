#!/usr/bin/env python3

"""
per_protein_protparam (+ metapredict)
-------------------------------------
Processes a directory of FASTA files (.faa) and, for each file:

1. Computes ProtParam metrics for each protein:
   - length
   - gravy
   - pI
   - percent_A..percent_Y

2. (Optional) Runs metapredict on the ENTIRE .faa file and computes
   six disorder summaries per protein:
   - dis_mean = mean disorder score
   - dis_frac_0_5 = % of residues with score >= 0.5
   - dis_max_run_0_5 = longest continuous stretch with score >= 0.5
   - dis_top10_mean = mean of the top 10% most disordered residues
   - dis_N_term30 = mean disorder of the first 30 amino acids
   - dis_C_term30 = mean disorder of the last 30 amino acids

3. Writes one TSV per .faa with rows = proteins and columns:
   seq_id, length, gravy, pI, percent_A..percent_Y,
   dis_mean, dis_frac_0_5, dis_max_run_0_5,
   dis_top10_mean, dis_N_term30, dis_C_term30
"""

# ---------------------------------------------------------------------
# 0. Imports
# ---------------------------------------------------------------------

import argparse                    # Read command-line arguments
import math
import multiprocessing as mp       # Parallelize processing across faa files
from pathlib import Path           # Robust path handling
from typing import Dict, List

import pandas as pd                # Build and save TSV files
from Bio import SeqIO              # Biopython: FASTA parsing
from Bio.SeqUtils.ProtParam import ProteinAnalysis  # pI, GRAVY, AA frequencies
from tqdm import tqdm              # Progress bar
import sys                         # Write warnings to stderr
import metapredict as meta         # Per-residue disorder prediction

STANDARD = list("ACDEFGHIKLMNPQRSTVWY")  # Standard amino acids
AMBIGUOUS = set("BJZXUO*-")              # Characters removed to avoid errors

# ---------------------------------------------------------------------
# 1. Sequence cleaning
# ---------------------------------------------------------------------

def clean_sequence(
    seq: str,
    strip_ambiguous: bool = True
) -> str:
    """
    Receives a sequence (string) and:
      - converts it to uppercase
      - removes spaces and line breaks
      - optionally removes ambiguous characters

    Returns the cleaned sequence.
    """
    s = str(seq).upper().replace(" ", "").replace("\n", "")

    if strip_ambiguous:
        s = "".join(ch for ch in s if ch not in AMBIGUOUS)

    return s


# ---------------------------------------------------------------------
# 2. ProtParam per protein
# ---------------------------------------------------------------------

def per_protein_metrics(
    faa_path: Path,
    strip_ambiguous: bool = True,
    keep_empty: bool = False
) -> pd.DataFrame:
    """
    Reads a .faa file and returns a DataFrame with ONE ROW per protein.

    Columns:
      seq_id, length, gravy, pI, percent_A..percent_Y
    """

    rows: List[Dict] = []

    # Iterate through all sequences in the .faa file
    for rec in SeqIO.parse(str(faa_path), "fasta"):

        raw = str(rec.seq)
        aa = clean_sequence(
            raw,
            strip_ambiguous=strip_ambiguous
        )

        if not aa:

            if keep_empty:
                # If the sequence becomes empty, fill with NaN and continue
                rows.append({
                    "seq_id": rec.id,
                    "length": 0,
                    "gravy": math.nan,
                    "pI": math.nan,
                    **{
                        f"percent_{aa}": math.nan
                        for aa in STANDARD
                    }
                })

            continue

        try:
            # Biopython ProtParam
            pa = ProteinAnalysis(aa)

            gravy = float(pa.gravy())                  # Mean hydropathy
            pI = float(pa.isoelectric_point())         # Isoelectric point

            # Use the newer property to avoid warnings
            perc = pa.amino_acids_percent

            row = {
                "seq_id": rec.id,
                "length": len(aa),
                "gravy": gravy,
                "pI": pI,
            }

            # Add the 20 amino acid frequencies
            for aa1 in STANDARD:
                row[f"percent_{aa1}"] = float(
                    perc.get(aa1, 0.0)
                )

            rows.append(row)

        except Exception:
            # If ProteinAnalysis fails, store NaNs
            rows.append({
                "seq_id": rec.id,
                "length": len(aa),
                "gravy": math.nan,
                "pI": math.nan,
                **{
                    f"percent_{a}": math.nan
                    for a in STANDARD
                }
            })

    # Convert list of dictionaries into a DataFrame
    df = pd.DataFrame(
        rows,
        columns=["seq_id", "length", "gravy", "pI"] +
                [f"percent_{a}" for a in STANDARD]
    )

    # Length is an inexpensive and highly useful metadata field
    # for downstream analyses and validation.

    return df


# ---------------------------------------------------------------------
# 3. METAPREDICT
# ---------------------------------------------------------------------

def normalize_header_to_seqid(header: str) -> str:
    """
    Metapredict returns the full FASTA header.

    Keep only the identifier portion, i.e. everything before '#'.
    """
    h = header.strip()

    if " # " in h:
        return h.split(" # ", 1)[0].strip()

    return h


def summarize_metapredict_scores(
    scores: List[float]
) -> Dict[str, float]:
    """
    Receives a list of disorder scores (one per residue, 0–1)
    and computes the six disorder metrics of interest.
    """

    L = len(scores)

    if L == 0:
        # No scores available: fill with NaN values
        return {
            "dis_mean": math.nan,
            "dis_frac_0_5": math.nan,
            "dis_max_run_0_5": math.nan,
            "dis_top10_mean": math.nan,
            "dis_N_term30": math.nan,
            "dis_C_term30": math.nan,
        }

    # Mean protein disorder
    dis_mean = sum(scores) / L

    # Fraction of clearly disordered residues
    thr = 0.5
    dis_frac_0_5 = sum(
        1 for s in scores if s >= thr
    ) / L

    # Longest run of disordered residues
    max_run = 0
    run = 0

    for s in scores:

        if s >= thr:
            run += 1

            if run > max_run:
                max_run = run

        else:
            run = 0

    # Mean of the top 10% most disordered residues
    # Avoids long proteins appearing weakly disordered overall
    k = max(1, int(0.10 * L))

    top_scores = sorted(
        scores,
        reverse=True
    )[:k]

    dis_top10_mean = (
        sum(top_scores) / len(top_scores)
    )

    # Mean disorder in N- and C-termini
    # (30 aa window, or shorter if protein length < 30)
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

        seq_id -> {
            dis_mean,
            dis_frac_0_5,
            ...
        }
    """

    pred_dict = meta.predict_disorder_fasta(
        str(faa_path)
    )

    out: Dict[str, Dict[str, float]] = {}

    for header, pair in pred_dict.items():

        seq_id = normalize_header_to_seqid(header)

        sequence_str, scores = pair

        # Sequence itself is not needed;
        # keep only summarized disorder statistics
        out[seq_id] = summarize_metapredict_scores(
            scores
        )

    return out


# ---------------------------------------------------------------------
# 4. Worker: process ONE .faa file
# ---------------------------------------------------------------------

def process_file(args_tuple):

    (
        faa_path,
        out_dir,
        strip_ambiguous,
        keep_empty,
        use_metapredict,
    ) = args_tuple

    # Convert strings to Path objects
    faa_path = Path(faa_path)
    out_dir = Path(out_dir)

    out_dir.mkdir(
        parents=True,
        exist_ok=True
    )

    # One TSV per .faa
    out_tsv = out_dir / f"{faa_path.stem}.protparam.tsv"

    # 1) ProtParam → DataFrame
    df = per_protein_metrics(
        faa_path,
        strip_ambiguous=strip_ambiguous,
        keep_empty=keep_empty
    )

    # 2) metapredict (optional)
    if use_metapredict:

        try:
            met_map = run_metapredict_on_faa(
                faa_path
            )

            # Add six disorder columns by mapping seq_id
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

    # 3) Save the final TSV
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


# ---------------------------------------------------------------------
# 5. main(): prepare jobs and launch Pool
# ---------------------------------------------------------------------

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
        "--no-strip-ambiguous",
        action="store_true",
        help="Do NOT remove ambiguous characters"
    )

    ap.add_argument(
        "--keep-empty",
        action="store_true",
        help="Write empty rows if sequences become empty"
    )

    ap.add_argument(
        "--no-metapredict",
        action="store_true",
        help="Do not run metapredict (ProtParam only)"
    )

    args = ap.parse_args()

    # Expand ~ in paths
    in_dir = Path(args.input_dir).expanduser()
    out_dir = Path(args.out_dir).expanduser()

    # List .faa files
    paths = sorted(
        in_dir.glob(args.glob)
    )

    if not paths:
        raise SystemExit(
            f"No files matched {in_dir}/{args.glob}"
        )

    # Cleaning and metapredict options
    strip_ambiguous = not args.no_strip_ambiguous
    use_metapredict = not args.no_metapredict

    # One job per .faa file
    jobs = [
        (
            p,
            out_dir,
            strip_ambiguous,
            args.keep_empty,
            use_metapredict
        )
        for p in paths
    ]

    if args.workers > 1:

        # Parallel processing using multiprocessing.Pool
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

    out_dir.mkdir(
        parents=True,
        exist_ok=True
    )

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
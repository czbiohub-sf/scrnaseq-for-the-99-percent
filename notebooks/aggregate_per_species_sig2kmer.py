
import argparse
import os
import glob
import re

from tqdm import tqdm
import pandas as pd
from joblib import Parallel, delayed
from sig2kmer import read_kmer_csv
from sig_utils import SKETCH_INFO_PATTERN
from sourmash.logging import error, notify, set_quiet
import sqlite3


def process_single_kmer_csv(csv):
    try:
        df = read_kmer_csv(csv)
    except pd.errors.EmptyDataError:
        #         print(f"Empty file: {csv}")
        return

    split = csv.split("/")
    test_species = re.findall("((test|train)-\w+)", csv)[0][0]
    species = test_species.split("-")[-1]
    df["species"] = species
    df["cell_id"] = os.path.basename(csv).split(".")[0]
    mol_or_alpha, moltype, ksize, style, value = re.findall(SKETCH_INFO_PATTERN, csv)[0]
    df["sketch_id"] = split[-4]
    df["moltype"] = moltype
    df["ksize"] = ksize
    df[style] = int(value)
    alignment_status = split[-2]
    df["alignment_status"] = alignment_status
    return df


def main():
    p = argparse.ArgumentParser()
    # base directory containing a 2--single-cell-kmers folder which contains sketch id directories with sig2kmer csvs
    p.add_argument("species_base_dir")
    p.add_argument(
        "--n-jobs",
        type=int,
        default=16,
        help="Number of processes to use",
    )
    p.add_argument(
        "--kmer-subdir",
        default="2--single-cell-kmers",
        type=str,
        help="Subdirectory containing csvs within each per-sketch id subdirectory",
    )
    p.add_argument(
        "--no-aligned-unaligned-subdir",
        action="store_true",
        help=(
            "If not set, looks for files in {species_base_dir}/{kmer_subdir}/{sketch_id}/csvs/aligned and "
            "{species_base_dir}/{kmer_subdir}/{sketch_id}/csvs/unaligned. Otherwise, "
            "looks in {species_base_dir}/{kmer_subdir}/{sketch_id}/csvs/"
        ),
    )
    args = p.parse_args()

    kmer_dir = os.path.join(args.species_base_dir, args.kmer_subdir)

    sketch_globber = os.path.join(
        kmer_dir,
        "alphabet-*ksize-*",
    )

    for sketch_dir in glob.glob(sketch_globber):
        notify(f"Reading hash2kmer csvs from {sketch_dir}")
        if args.no_aligned_unaligned_subdir:
            csv_globber = os.path.join(
                sketch_dir,
                "csvs",
                "*.csv",
            )
        else:
            csv_globber = os.path.join(
                sketch_dir,
                "csvs",
                "*",  # "aligned" or "unaligned" directory
                "*.csv",
            )
        total = sum(1 for _ in glob.iglob(csv_globber))

        dfs = Parallel(n_jobs=args.n_jobs)(
            delayed(process_single_kmer_csv)(csv) for csv in glob.iglob(csv_globber)
        )
        try:
            kmers = pd.concat(dfs)
        except ValueError:
            # No objects to contatenate, continue
            notify("No hash2kmer files found, skipping")
            continue
        parquet = os.path.join(sketch_dir, "hash2kmer.parquet")
        notify(f"Writing {parquet}")
        kmers.to_parquet(parquet)


if __name__ == "__main__":
    main()

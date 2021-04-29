
import argparse
import glob
import os

import pandas as pd
import scanpy as sc
from joblib import Parallel, delayed
from IPython.display import display
from tqdm import tqdm

SHARED_CELLTYPES = [
    "Alveolar Epithelial Type 2",
    "B cell",
    "Capillary",
    "Dendritic",
    "Fibroblast",
    "Macrophage",
    "Monocyte",
    "Natural Killer T cell",
    "Smooth Muscle and Myofibroblast",
    "T cell",
]


def describe(df, random=False):
    print(df.shape)
    print("--- First 5 entries ---")
    display(df.head())
    if random:
        print("--- Random subset ---")
        display(df.sample(5))


def process_hash2kmer(parquet, adata_shared, celltype_col):
    hash2kmer = pd.read_parquet(parquet)
    describe(hash2kmer)

    hash2kmer_with_celltypes = hash2kmer.join(
        adata_shared.obs[celltype_col], on="cell_id"
    )

    hash2kmer_celltype_unique_hashvals = hash2kmer_with_celltypes.drop_duplicates(
        [
            "kmer_in_sequence",
            "kmer_in_alphabet",
            "hashval",
            "gene_name",
            "alignment_status",
            "broad_group",
        ]
    )
    describe(hash2kmer_celltype_unique_hashvals)

    parquet_out = parquet.replace(".parquet", "__unique_kmers_per_celltype.parquet")
    hash2kmer_celltype_unique_hashvals.to_parquet(parquet_out)

    # Show number of aligned/unaligned k-mers per celltype
    per_celltype_alignment_status_kmers = hash2kmer_celltype_unique_hashvals.groupby(
        celltype_col, observed=True
    ).alignment_status.value_counts()
    print(per_celltype_alignment_status_kmers)


def main():
    p = argparse.ArgumentParser()
    # base directory containing a 2--single-cell-kmers folder which contains sketch id directories with sig2kmer csvs
    p.add_argument("species_base_dir")
    p.add_argument(
        "--kmer-subdir",
        default="2--single-cell-kmers",
        type=str,
        help="Subdirectory containing csvs within each per-sketch id subdirectory",
    )
    p.add_argument(
        "--h5ad",
        default="/home/olga/data_sm/immune-evolution/h5ads/human-lemur-mouse-bat/human-lemur-mouse-bat__lung_only.h5ad",
        help=("Location of the AnnData h5ad object of single-cell data"),
    )
    p.add_argument(
        "--n-jobs",
        default=3,
        type=int,
        help=(
            "Number of jobs to do in parallel. By default, 3 for the 3 molecule types (DNA, protein, Dayhoff)"
        ),
    )
    p.add_argument(
        "--celltype-col",
        default="broad_group",
        help=(
            "Column name endcoding the cell type in the h5ad AnnData object, i.e. an adata.obs column"
        ),
    )

    args = p.parse_args()

    adata = sc.read(args.h5ad)
    adata.obs = adata.obs.reset_index().set_index("cell_id")

    adata_shared = adata[adata.obs[args.celltype_col].isin(SHARED_CELLTYPES)]

    parquets = glob.iglob(
        os.path.join(
            args.species_base_dir,
            args.kmer_subdir,
            "*",  # This is the sketch_id, e.g. alphabet-DNA__ksize-21__scaled-10
            "hash2kmer.parquet",
        )
    )

    if args.n_jobs > 1:
        Parallel(n_jobs=args.n_jobs)(
            delayed(process_hash2kmer)(parquet, adata_shared, args.celltype_col)
            for parquet in parquets
        )
    else:
        for parquet in tqdm(parquets):
            print("hash2kmer parquet:", parquet)
            process_hash2kmer(parquet, adata_shared, args.celltype_col)


if __name__ == "__main__":
    main()

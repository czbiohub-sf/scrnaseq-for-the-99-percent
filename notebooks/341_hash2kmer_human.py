#!/usr/bin/env python
# coding: utf-8


import glob
import os
import re

import pandas as pd
import scanpy as sc
import sig_utils
import sourmash
from joblib import Parallel, delayed
from tqdm import tqdm

from pandas_utils import describe


def remove_rogue_tqdm():
    import tqdm

    try:
        tqdm._instances.clear()
    except AttributeError:
        pass


output_dir = "/home/olga/data_sm/immune-evolution/kmer-signatures/1--train-mouse/5-celltype-kmers"


def describe(df, random=False):
    print(df.shape)
    print("--- First 5 entries ---")
    display(df.head())
    if random:
        print("--- Random subset ---")
        display(df.sample(5))


SHARED_CELLTYPES = [
    "Capillary",
    "Alveolar Epithelial Type 2",
    "B cell",
    "T cell",
    "Natural Killer T cell",
    "Macrophage",
    "Monocyte",
    "Dendritic",
    "Fibroblast",
    "Smooth Muscle and Myofibroblast",
]


sanitized_to_shared = {sig_utils.sanitize(x): x for x in SHARED_CELLTYPES}
sanitized_to_shared


h5ad = os.path.join(
    "/home/olga/data_sm/immune-evolution/data-objects/human-lemur-mouse-bat/",
    "human-lemur-mouse-bat__lung_only.h5ad",
)
adata = sc.read(h5ad)
adata.obs = adata.obs.set_index("cell_id")


mouse_kmermaid_base = "/home/olga/data_sm/immune-evolution/pipeline-results/mouse/kmermaid/lung--mouse--remove-ribo"


mouse_translate_folder = os.path.join(mouse_kmermaid_base, "translate")


train_3_merged_celltype_remove_common = "/home/olga/data_sm/immune-evolution/kmer-signatures/1--train-mouse/3--merged-celltype-remove-common-kmers"

aligned_unaligned = "aligned", "unaligned"

for folder in glob.glob(
    os.path.join(train_3_merged_celltype_remove_common, "alphabet-*")
):

    sketch_id = os.path.basename(folder)

    sketch_id_output_dir = os.path.join(
        output_dir,
        sketch_id,
    )
    seqout_dir = os.path.join(sketch_id_output_dir, "fastas")
    for alignment_status in aligned_unaligned:
        os.makedirs(os.path.join(seqout_dir, alignment_status))

    sketch_info = re.findall(sig_utils.SKETCH_INFO_PATTERN, sketch_id)[0]
    moltype = sketch_info[1]
    ksize = int(sketch_info[2])
    scaled = int(sketch_info[4])
    print(sketch_id, sketch_info)

    if moltype == "DNA":
        fasta_moltype = "nucleotides"
    else:
        fasta_moltype = "peptides"

    glob_path = os.path.join(
        mouse_translate_folder, f"*__coding_reads_{fasta_moltype}.fasta"
    )

    total = sum(1 for _ in glob.iglob(glob_path))
    print(f"total number of fastas: {total}")

    sigfiles = glob.glob(os.path.join(folder, "*.sig"))

    for sigfile in sigfiles:
        query_hashvals = set([])
        sanitized = os.path.basename(sigfile).split(".sig")[0]
        celltype = sanitized_to_shared[sanitized]

        adata_subset = adata[adata.obs.broad_group == celltype]

        sigobjs = sourmash.load_file_as_signatures(
            sigfile, ksize=ksize, select_moltype=moltype
        )
        for sigobj in sigobjs:
            query_hashvals.update(sigobj.minhash.hashes.keys())

        dfs = Parallel(n_jobs=96)(
            delayed(sig_utils.get_hashvals_in_singlecell_fasta)(
                singlecell_fasta,
                cell_ontology_classes=adata_subset.obs["broad_group"],
                seqout_dir=seqout_dir,
                query_hashvals=query_hashvals,
                double_aligned=True,
                ksize=ksize,
                moltype=moltype,
            )
            for singlecell_fasta in tqdm(glob.iglob(glob_path), total=total)
        )
        celltype_kmers = pd.concat(dfs)
        celltype_kmers["celltype"] = celltype
        celltype_kmers["celltype_sanitized"] = sanitized
        describe(celltype_kmers)
        parquet = os.path.join(
            sketch_id_output_dir, f"hashvals-to-kmers__{sanitized}.parquet"
        )
        celltype_kmers.to_parquet(parquet)




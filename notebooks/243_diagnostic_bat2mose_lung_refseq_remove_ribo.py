#!/usr/bin/env python
# coding: utf-8

# # Imports

# In[1]:


import os

from joblib import Parallel, delayed
import pandas as pd
import sig_utils
from tqdm import tqdm
import sig2kmer


# ## Constants (ksizes, scaled)


STYLE = "scaled"
SCALED = 10
ALPHA = "dayhoff"
KSIZE = 51

SKETCH_ID = sig_utils.make_sketch_id(
    alpha=ALPHA, ksize=KSIZE, style=STYLE, value=SCALED
)
SKETCH_ID


# ## Base dir for output


base_dir = (
    "/home/olga/data_sm/immune-evolution/analyses/2021-02-01__search_bat_in_mouse"
)


# ## Def describe


def describe(df, random=False):
    print(df.shape)
    print("--- First 5 entries ---")
    display(df.head())
    if random:
        print("--- Random subset ---")
        display(df.sample(5))


# ## mouse and human sig paths

# ### Mouse sig paths


mouse_kmermaid_base = "/home/olga/data_sm/immune-evolution/pipeline-results/mouse/kmermaid/lung--mouse--remove-ribo"

bat_kmermaid_base = (
    "/home/olga/data_sm/immune-evolution/pipeline-results/bat/kmermaid/ksize_8"
)

mouse_sig_output_folder = os.path.join(
    bat_kmermaid_base,
    "search_bat_in_mouse_no_ribosome_no_containment_no_dissociation__per_celltype_remove_common_hashes__remove_refseq_ribosomal"
    #     "sketches_peptide_handmade_mouse_merged_remove_ribosomal_dissociation_all_mouse_cells",
)

mouse_per_celltype_sig_folder = os.path.join(
    bat_kmermaid_base,
    "sketches_peptide_handmade_mouse_merged_remove_ribosomal_dissociation_all_mouse_cells_per_celltype_remove_shared_hashes__refseq_remove_ribo",
)


# ### Bat sig paths


bat_sig_base = os.path.join(
    bat_kmermaid_base,
    "sketches_peptide_handmade_merged_downsampled"
    #     "sketches_peptide_handmade_merged_remove_ribosomal_dissociation_allcells",
)


bat_kmermaid_v2_base = "/mnt/ibm_sm/olga/immune-evolution/pipeline-results/bat/kmermaid/redo_singlecell_fastas"
bat_translate_dir = os.path.join(
    bat_kmermaid_v2_base,
    "translate",
)

# Make sure these paths exist and have stuff
get_ipython().system(" ls -lha $bat_translate_dir | head")


# ## Read parquet file of top hits

# In[13]:


parquet_fn = "sourmash_search_bat_in_mouse_remove_ribosomal_dissociation__predicted_cells_top_hit.parquet"
parquet = os.path.join(base_dir, parquet_fn)
predicted_cells_top_hit = pd.read_parquet(parquet)
predicted_cells_top_hit = predicted_cells_top_hit.query(
    "(ksize == @KSIZE) and (alphabet == @ALPHA)"
)
describe(predicted_cells_top_hit)

# ## Actually run diagnostic kmers

# In[14]:


base_dir


# In[15]:


seqout_dir = os.path.join(base_dir, "diagnostic_kmer_seqs_bat")
get_ipython().system(" rm -rf $seqout_dir")
get_ipython().system(" mkdir -p $seqout_dir/aligned")
get_ipython().system(" mkdir -p $seqout_dir/unaligned")
print(seqout_dir)


# In[16]:


dfs = []

alignment_statuses = "aligned", "unaligned"

# for i, row in tqdm(predicted_cells_top_hit.iterrows()):


def get_diagnostic_kmers_per_row(
    cell_id,
    predicted_celltype_name,
    ground_truth_celltype_name,
    seqout_dir,
    query_sig_base,
    query_translate_dir,
    db_sig_folder,
    double_aligned=False,
    channel_suffix="_possorted_genome_bam",
    add_ksize_to_sig_path=True,
):

    pred_sanitized = sig_utils.sanitize(predicted_celltype_name)
    true_sanitized = sig_utils.sanitize(ground_truth_celltype_name)
    seqout_template = os.path.join(
        seqout_dir,
        "{alignment_status}",
        sig_utils.sanitize(predicted_celltype_name),
        f"pred-{pred_sanitized}---{cell_id}---true-{true_sanitized}.fasta",
    )

    df = sig2kmer.get_diagnostic_kmers_for_cell(
        cell_id,
        cell_sig_base=query_sig_base,
        cell_fasta_dir=query_translate_dir,
        predicted_celltype_name=predicted_celltype_name,
        celltype_sig_base=db_sig_folder,
        alphabet=ALPHA,
        ksize=KSIZE,
        scaled=SCALED,
        double_aligned=double_aligned,
        verbose=False,
        seqout_template=seqout_template,
        channel_suffix=channel_suffix,
        add_ksize_to_sig_path=True,
    )
    return df


def parallel_get_all_diagnostic_kmers(
    predicted_cells_top_hit,
    query_sig_base,
    query_translate_dir,
    predicted_celltype_col,
    groundtruth_celltype_col,
):
    dfs = Parallel(n_jobs=96)(
        delayed(get_diagnostic_kmers_per_row)(
            cell_id=row.name,
            predicted_celltype_name=row[predicted_celltype_col],
            ground_truth_celltype_name=row[groundtruth_celltype_col],
            seqout_dir=seqout_dir,
            query_sig_base=query_sig_base,
            query_translate_dir=query_translate_dir,
            mouse_per_celltype_sig_folder=mouse_per_celltype_sig_folder,
            double_aligned=True,
            channel_suffix="_possorted_genome_bam",
        )
        for (i, row) in tqdm(
            predicted_cells_top_hit.iterrows(), total=len(predicted_cells_top_hit)
        )
    )
    diagnostic_kmers = pd.concat(dfs)
    describe(diagnostic_kmers)
    return diagnostic_kmers


parallel_get_all_diagnostic_kmers()

# In[20]:


diagnostic_kmers["alphabet"] = ALPHA
diagnostic_kmers["ksize"] = KSIZE


# In[21]:


# ## Write diagnostic kmers to file

# In[22]:


parquet = os.path.join(
    base_dir,
    f"diagnostic_kmers_bat_with_gene_names__alpha-{ALPHA}__ksize-{KSIZE}.parquet",
)
diagnostic_kmers.to_parquet(parquet)
# In[20]:


diagnostic_kmers.groupby(
    ["predicted_celltype_name", "alignment_status"]
).hashval.nunique()


# In[21]:


n_kmers_per_celltype = diagnostic_kmers.groupby(
    ["predicted_celltype_name", "alignment_status"]
).hashval.nunique()
n_kmers_per_celltype


# In[22]:


percent_kmers_per_celltype = n_kmers_per_celltype.groupby(level=0).apply(
    lambda x: 100 * x / x.sum()
)
percent_kmers_per_celltype.name = "percent_cells"
percent_kmers_per_celltype


# In[23]:


csv = os.path.join(base_dir, "diagnostic_kmers_bat__percent_cells_per_celltype.csv")
percent_kmers_per_celltype.to_csv(csv)


# In[ ]:


# In[24]:


diagnostic_kmer_genes = diagnostic_kmers.groupby(
    ["predicted_celltype_name", "alignment_status", "gene_name"]
).abundances.sum()
diagnostic_kmer_genes = diagnostic_kmer_genes.groupby(level=[0, 1, 2]).nlargest(5)
# diagnostic_kmer_genes = diagnostic_kmer_genes[diagnostic_kmer_genes > 10]
diagnostic_kmer_genes


# In[25]:


top5_genes_per_celltype = diagnostic_kmer_genes.groupby(
    level=[0, 1], group_keys=False
).apply(lambda x: x.sort_values(ascending=False).head(5))
top5_genes_per_celltype


# In[26]:


top5_genes_per_celltype_df = top5_genes_per_celltype.reset_index()
loc_genes = top5_genes_per_celltype_df.gene_name[
    top5_genes_per_celltype_df.gene_name.str.startswith("LOC")
].unique()
loc_genes


# In[27]:


diagnostic_kmers_loc_genes = diagnostic_kmers.query("gene_name in @loc_genes")
describe(diagnostic_kmers_loc_genes)


# In[28]:


diagnostic_kmers_loc_genes.groupby("gene_name").hashval.nunique()


# In[29]:


loc_gene_hashes = diagnostic_kmers_loc_genes.groupby("gene_name").hashval.unique()
loc_gene_hashes


# In[42]:


loc_gene_hashes["LOC109443676"]


# In[ ]:


set(diagnostic_kmers_loc_genes.hashval)


# In[ ]:

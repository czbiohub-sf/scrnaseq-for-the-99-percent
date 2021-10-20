
import os

FIGURE_FOLDER = os.path.join('..', 'figures')
DATA_FOLDER = os.path.join('..', 'kmer-homology-data')


sig_outdir_base = "/home/olga/data_lg/data_sm_copy/immune-evolution/kmer-signatures"

top_hit_suffix = os.path.join(
    "4--aggregated-results",
    "sourmash-search-results--top-hit.parquet",
)

self2self_parquet = os.path.join(
    sig_outdir_base, "0--mouse2mouse", "0--self2self-bootstrapped", top_hit_suffix
)

mouse2mouse_parquet = os.path.join(
    sig_outdir_base, "0--mouse2mouse", "1--mouse2mouse", top_hit_suffix
)

lemur_parquet = os.path.join(sig_outdir_base, "4--test-lemur", top_hit_suffix)

bat_parquet = os.path.join(sig_outdir_base, "3--test-bat", top_hit_suffix)
human_parquet = os.path.join(sig_outdir_base, "2--test-human", top_hit_suffix)

top_hit_paths = {
#     "self": self2self_parquet,
    "mouse": mouse2mouse_parquet,
    "lemur": lemur_parquet,
    "bat": bat_parquet,
    "human": human_parquet,
}

H5AD = os.path.join(DATA_FOLDER, 'h5ads',
 'human-lemur-mouse-bat',
 'human-lemur-mouse-bat__lung_only.h5ad')
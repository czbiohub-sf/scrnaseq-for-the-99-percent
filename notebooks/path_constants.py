import os

ORPHEUM_BENCHMARKING_FOLDER = "/Users/olgabot/Downloads/00-orpheum-benchmarking"
ORPHEUM_EXTERNAL_DATA = os.path.join(ORPHEUM_BENCHMARKING_FOLDER, "00-external-data")

QFO_EUKARYOTA_FOLDER = os.path.join(
    ORPHEUM_EXTERNAL_DATA, "quest-for-orthologs", "2019", "Eukaryota"
)


ORPHEUM_GROUND_TRUTH_FOLDER = os.path.join(
    ORPHEUM_BENCHMARKING_FOLDER, "02-results", "00-ground-truth-protein-coding-frames"
)

SIMULATED_RNASEQ_FOLDER = os.path.join(
    ORPHEUM_BENCHMARKING_FOLDER, "01-simulated-rnaseq"
)

MAMMALIA_BUSCO_SUBSET_FOLDER = os.path.join(
    ORPHEUM_BENCHMARKING_FOLDER, "mammalia_busco_subsets"
)


FIGURE_FOLDER = os.path.join("..", "figures")
DATA_FOLDER = os.path.join("..", "kmer-homology-data")


FIGURE_FOLDER = os.path.join("..", "figures")
DATA_FOLDER = os.path.join("..", "kmer-homology-data")

# Raw data = straight from the sequencer, or downloaded from a database
RAWDATA_FOLDER = os.path.join(DATA_FOLDER, "00--rawdata")

### --- Figure 2, Supplementary Figure 1 Paths --- ###

# Quest for Orthologs 2019 folder
QFO_FOLDER = os.path.join(RAWDATA_FOLDER, "quest-for-orthologs", "2019")

# Quest for Orthologs 2019 Eukaryota data
QFO_EUKARYOTA_FOLDER = os.path.join(
    RAWDATA_FOLDER, "quest-for-orthologs", "2019", "Eukaryota"
)

# Ortholog database (OrthoDB) folder
ORTHODB_FOLDER = os.path.join(RAWDATA_FOLDER, "orthodb", "v10.1")

# BUSCO Mammalia folder
BUSCO_MAMMALIA_FOLDER = os.path.join(RAWDATA_FOLDER, "busco", "mammalia_odb10", "info")

# Processed data = data that has been maniuplated in some way from the raw data
PROCESSED_DATA_FOLDER = os.path.join(DATA_FOLDER, "01--processed-data")


# Data related to orpheum benchmarking
ORPHEUM_BENCHMARKING_FOLDER = os.path.join(
    PROCESSED_DATA_FOLDER, "orpheum-benchmarking"
)


# Human simulated reads -- output from Polyster simulated reads
SIMULATED_READS_FOLDER = os.path.join(ORPHEUM_BENCHMARKING_FOLDER, "simulated", "human")

# Folder of output of Orpheum translate on human reads
ORPHEUM_PIPELINE_RESULTS_FOLDER = os.path.join(
    ORPHEUM_BENCHMARKING_FOLDER, "orpheum-translate-on-simulated-human-data-results"
)

### --- Figure 1, Supplementary Figure 2 Paths --- ###

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

H5AD = os.path.join(
    DATA_FOLDER,
    "h5ads",
    "human-lemur-mouse-bat",
    "human-lemur-mouse-bat__lung_only.h5ad",
)

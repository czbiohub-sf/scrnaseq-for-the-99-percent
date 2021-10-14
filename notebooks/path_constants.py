import os

FIGURE_FOLDER = os.path.join("..", "figures")
DATA_FOLDER = os.path.join("..", "kmer-homology-data")

# Raw data = straight from the sequencer, or downloaded from a database
RAWDATA_FOLDER = os.path.join(DATA_FOLDER, "00--rawdata")

# Processed data = data that has been maniuplated in some way from the raw data
PROCESSED_DATA_FOLDER = os.path.join(DATA_FOLDER, "01--processed-data")

# Data related to orpheum benchmarking
ORPHEUM_BENCHMARKING_FOLDER = os.path.join(
    PROCESSED_DATA_FOLDER, "orpheum-benchmarking"
)
MAMMALIA_BUSCO_SUBSET_FOLDER = os.path.join(
    ORPHEUM_BENCHMARKING_FOLDER, "mammalia_busco_subsets"
)

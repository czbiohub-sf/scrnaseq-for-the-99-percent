import re
import json
import argparse
import pandas as pd
from pandas import json_normalize

#ncRNA = glob.glob("/home/ntpierce/2020-simulate-rnaseq/output_simreads/sencha/translate/*ncRNA*json")
#qfo_dna = glob.glob("/home/ntpierce/2020-simulate-rnaseq/output_simreads/sencha/translate/*qfo_dna*json")
#species_metadata = pd.read_csv("/home/ntpierce/2020-simulate-rnaseq/species_metadata.csv")

def build_summary_df(filename_list):
    all_info=[]
    for fi in filename_list:
        with open(fi) as f:
            all_info.append(json_normalize(json.load(f)))
    return all_info

def find_origin(row):
    if "qfo" in row["input_files"][0]:
        row["origin"] = "coding"
        row["true_positives"] = row["categorization_counts.Coding"]
        row["true_negatives"] = 0
        row["false_positives"] = 0
        row["false_negatives"] = row["total_reads"] - row["categorization_counts.Coding"]
        input_file = row["input_files"][0]
        row["species"] = re.search("(\w*)_qfo", input_file).groups()[0].replace("_", " ")
    elif "ncRNA" in row["input_files"][0]:
        row["origin"] = "noncoding"
        row["true_positives"] = 0
        row["true_negatives"] = row["total_reads"] - row["categorization_counts.Coding"]
        row["false_positives"] = row["categorization_counts.Coding"]
        row["false_negatives"]= 0
        input_file = row["input_files"][0]
        row["species"] = re.search("(\w*)_ncRNA", input_file).groups()[0].replace("_", " ")
        #row["species"] = row["input_files"].str.extract("(\w*)_ncRNA")[0]
    if "Hsapiens" in row["peptide_bloom_filter"]:
        row["peptide_reference"] = "Homo sapiens QfO"
    elif "sprot" in row["peptide_bloom_filter"]:
        row["peptide_reference"] = "Swiss Prot"
    return row


def summarize_json(args):
    input_files = args.input_files
    out_csv = args.output_csv
    species_metadata_file = args.species_metadata

    # build data frame from input json files
    summaryInfo = build_summary_df(input_files)
    summaryDF = pd.concat(summaryInfo).fillna(0)
    # calculate and add some metrics: total_reads, origin (coding/noncoding), TP,FP,TN,FN. bespoke to my filenaming for now
    summaryDF["total_reads"] = summaryDF.groupby((summaryDF.columns.str.split(".").str[0].str.contains("categorization_counts")),axis=1).sum()[True]
    summaryDF = summaryDF.apply(find_origin, axis=1)
    # read in species metadata, add divergence times to dataframe
    species_metadata = pd.read_csv(species_metadata_file)
    summaryDF=summaryDF.merge(species_metadata[["scientific_name", "divergence_from_human_mya"]], left_on=["species"], right_on=["scientific_name"])

    # groupby the relevant columns and sum true pos, true neg, false pos, false negs
    # this is important bc we have two files for each species: one of known coding genes, one of known noncoding sequences
    # thus each file will set two of these metrics (see find_origin), and we need to sum over both files per species to get all four.
    metricDF = summaryDF.groupby(["species", "divergence_from_human_mya", "peptide_alphabet", "peptide_reference", "peptide_ksize", "peptide_bloom_filter","jaccard_threshold"]).sum()[["true_positives", "true_negatives", "false_positives", "false_negatives"]]

    # now, with all tp/tn/fp/fn, calculate precision, recall, and f1_score
    metricDF["precision"] = metricDF["true_positives"]/(metricDF["true_positives"] + metricDF["false_positives"])
    metricDF["recall"] = metricDF["true_positives"]/(metricDF["true_positives"] + metricDF["false_negatives"])
    metricDF["F1_score"] = 2 * ((metricDF["precision"] * metricDF["recall"]) / (metricDF["precision"] + metricDF["recall"]))
    # reset the index
    metricDF = metricDF.reset_index()
    metricDF.to_csv(out_csv, index=False)


if __name__ == "__main__":
    p = argparse.ArgumentParser()
    p.add_argument("input_files", nargs='+')
    p.add_argument("--species-metadata", default="species_metadata.csv", help="provide metadata file with species divergence times")
    p.add_argument("--output-csv", help="name for output csv")
    args = p.parse_args()
    sys.exit(summarize_json(args))


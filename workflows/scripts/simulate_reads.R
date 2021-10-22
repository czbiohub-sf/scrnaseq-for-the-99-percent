# tiny sim script: given fasta, simulate some reads using polyester

# run this script from within polyester-env.yml
log <- file(snakemake@log[[1]], open="wt")
sink(log)
sink(log, type="message")

library(polyester)
library(Biostrings)
library(rtracklayer)
library(BSgenome)

# Use fasta of interest to simulate reads

subset_fasta <- snakemake@input[[1]]
output_dir <- snakemake@params[["output_dir"]]
print(output_dir)
num_reps <- as.integer(snakemake@params[["num_reps"]]) # 5
print(num_reps)
read_length <- as.integer(snakemake@params[["read_length"]]) # 150
print(read_length)
paired <- snakemake@params[["simulate_paired"]] # true/false
print(paired)

reads_per_tx <- as.integer(snakemake@params[["num_reads_per_transcript"]])
print(reads_per_tx)

# enable either coverage or solid number
#coverage <- snakemake@params[["coverage"]]
#coverage = round(20 * width(fastaStrSet) / 100) # 20x coverage

# build count matrix
fastaStrSet <- readDNAStringSet(subset_fasta)
head(fastaStrSet)
countmat = matrix(reads_per_tx, nrow=length(readDNAStringSet(subset_fasta)), ncol=num_reps)
print(countmat)
# simulate reads
simulate_experiment_countmat(fasta=subset_fasta, readmat=countmat, reportCoverage=TRUE, readlen=read_length, paired=paired, outdir=output_dir, gzip=TRUE)

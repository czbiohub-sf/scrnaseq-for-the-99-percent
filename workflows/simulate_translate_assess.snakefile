"""
Author: N Tessa Pierce, UC Davis Lab for Data Intensive Biology
Run: snakemake -s simulate_translate_assess.snakefile --use-conda --configfiles config/QfO_vertebrates.yml
"""
# simulate RNAseq reads from transcripts in a fasta file, translate them using sencha translate

import os

envs_dir="envs"
outdir = config.get("output_directory","output_simreads")
logs_dir = os.path.join(outdir, "logs")

datadir_default = os.path.join(outdir, "data")
datadir = config.get(outdir, datadir_default)
simdir = os.path.join(outdir, "simulate_reads")
index_dir = os.path.join(outdir, "sencha_index")
translate_dir = os.path.join(outdir, "sencha_translate")

# use config to specify files
runname = config["run_name"]
samplefiles = config["samples"]
# download and simulate these types per species
default_sampletypes = ["ncRNA", "cds"]
sample_types = config.get("sample_types", default_sampletypes)

# if using different types of reference files (e.g. noncoding rna vs coding seqs), build a refname: download link dictionary
sampleLinks={}
for sample in samplefiles:
    for sampletype, link in samplefiles[sample].items():
        if sampletype in sample_types:
            samplename = f"{sample}_{sampletype}"
            sampleLinks[samplename] = link # build name : download link dictionary

## polyester creates files called "sample_01.fasta.gz", sample_02.fasta.gz", etc
def generate_replist(number_replicates):
    replist = list(range(number_replicates))
    reps=[]
    for num in replist:
        num+=1 # start at 1
        if num < 10:
            num = "0" + str(num)
        reps.append(str(num))
    return reps

# get polyester parameters from config; use replicates to build rep sample numbers
sim_config = config.get("simulation_params")
replist = generate_replist(sim_config.get("num_replicates", 5))


### build sencha translate outputs (need to loop through refs + tablesize, build targs based on alpha, ksize, jaccardthresh
reference_info = config["sencha_params"].get("peptide_references")
refnames = []
default_tablesize = "1e8"

# include tablesize in reference name
for ref, info in reference_info.items():
    tablesize = info.get("tablesize", default_tablesize)
    refnames.append(f"{ref}_t{tablesize}")

# build output files based on config
translate_json_targets, ref_targets = [],[]
#output_extensions = ["codingpep.fa", "noncoding.fa", "lowcomplexnucl.fa", "csv", "json"]

alphabets_to_run = config["sencha_params"]["alphabet"] # dictionary of alphabets to run, corresponding ksizes
# build final file targets based on config info
for alpha, alpha_info in alphabets_to_run.items():
    # sencha index targets
    ref_targets+=expand(os.path.join(index_dir, "ref{refname}_{alphabet}_k{k}.index"), refname=refnames, alphabet=alpha, k=alpha_info["ksizes"])
    # sencha translate targets
    translate_json_targets+=expand(os.path.join(translate_dir, "{sample}_{rep}.{alphabet}_k{k}_ref{refname}_jacc{thresh}.json"), sample=sampleLinks.keys(), alphabet=alpha, k=alpha_info["ksizes"], refname=refnames, thresh=alpha_info["jaccard_threshold"], rep=replist)
    

rule all:
    input: 
         expand(os.path.join(translate_dir, "{run_name}_summary.csv"), run_name=runname)
#        expand(os.path.join(simdir, "{run_name}_samples.csv"), run_name = runname), # this just makes csv 

rule download_fasta:
    output: os.path.join(datadir, "{samplename}.fa.gz")
    params: 
        download_link=lambda w: sampleLinks[w.samplename]
    log: os.path.join(logs_dir, "get_data", "{samplename}.log")
    threads: 1
    resources:
      mem_mb=1000, #1GB
      runtime=60 #minutes
    shell:
        """
        curl -L {params.download_link} -o {output} 2> {log}
        """

rule polyester_simreads_full:
    #input: lambda w: os.path.join(datadir, samplefiles[w.samplename]["fasta"])
    #input: lambda w: os.path.join(datadir,"{samplename}.fa.gz") 
    input: rules.download_fasta.output
    output: expand(os.path.join(simdir, "{{samplename}}", "sample_{rep}.fasta.gz"), rep = replist)
    params:
        output_dir = lambda w: os.path.abspath(os.path.join(simdir, w.samplename)),
        num_reps = sim_config.get("num_replicates", 5),
        read_length = sim_config.get("read_length", 150),
        simulate_paired = sim_config.get("paired", False),
        num_reads_per_transcript=sim_config.get("num_reads_per_transcript", 1000),
    log: os.path.join(logs_dir, "{samplename}.simreads.log")
    benchmark: os.path.join(logs_dir, "{samplename}.simreads.benchmark")
    threads: 1
    resources:
      mem_mb=16000, #16GB
      runtime=1000
    wildcard_constraints:
        samplename="\w+",
        gene_name="\w+"
    conda: "envs/polyester-env.yml"
    script: "scripts/simulate_reads.R"


rule seqtk_fasta_to_fastq_full:
    input: os.path.join(simdir, "{samplename}", "sample_{rep}.fasta.gz")
    output: os.path.join(simdir, "{samplename}", "{samplename}_{rep}.fq.gz")
    log: os.path.join(logs_dir, "seqtk", "{samplename}_{rep}.seqtk.log")
    benchmark: os.path.join(logs_dir, "seqtk", "{samplename}_{rep}.seqtk.benchmark")
    threads: 1
    resources:
      mem_mb=4000, #4GB
      runtime=60 #minutes
    wildcard_constraints:
        samplename="\w+",
        gene_name="\w+",
        rep="\d+"
    conda: "envs/seqtk-env.yml"
    shell:
        """
        seqtk seq -F 'I' {input} | gzip -9 > {output} 2> {log}
        """

#rule write_samples_csv:
#    input: expand(os.path.join(simdir, "{samplename}", "{samplename}_{rep}.fq.gz"), ref=sampleLinks.keys(), rep=replist),
#    output: os.path.join(simdir, "{run_name}_samples.csv")
#    threads: 1
#    resources:
#      mem_mb=1000,
#      runtime=15
#    wildcard_constraints:
#        samplename="\w+",
#        gene_name="\w+",
#        rep="\d+",
#        numreads="\d+"
#    run:
#        with open(str(output), "w") as csv:
#            csv.write("sample_id" + "," + "read1"+ "\n")
#            for f in input:
#                sample = os.path.basename(os.path.dirname(f))
#                replicate = f.rsplit(".fq.gz")[0].rsplit("_", 1)[1] # _01.fq.gz
#                csv.write(f"{samplename}_{replicate}" + "," + os.path.abspath(f) + "\n")

# abbreviate "hydrophobic-polar" as hp in filenames. Need full alpha name for sencha code
alpha_abbreviations = {"hp": "hydrophobic-polar", "protein": "protein", "dayhoff": "dayhoff"}

rule sencha_index:
    input: lambda w: reference_info[w.ref]["path"] # get file path from reference name, config
    output: os.path.join(index_dir, "ref{ref}_t{tablesize}_{alphabet}_k{ksize}.index"),
    log: os.path.join(logs_dir, "sencha_index", "ref{ref}_t{tablesize}_{alphabet}_k{ksize}.index.log")
    benchmark: os.path.join(logs_dir, "sencha_index", "ref{ref}_t{tablesize}_{alphabet}_k{ksize}.index.benchmark")
    params:
        alphabet=lambda w: alpha_abbreviations[w.alphabet],
    threads: 1
    resources:
        mem_mb=lambda w:  reference_info[w.ref].get("index_memory", 20000),
        runtime=1000,
    wildcard_constraints:
        ref="\w+",
        alphabet="\w+",
        ksize="\d+",
    conda: os.path.join(envs_dir, "sencha-env.yml")
    shell:
        """
        sencha index --alphabet {params.alphabet} --peptide-ksize {wildcards.ksize} {input} \
        --tablesize {wildcards.tablesize} --save-as {output} 2> {log}
        """

rule sencha_translate:
    input:
        #fastq=os.path.join(simdir, "{samplename}", "{samplename}_{rep}.fq.gz")
        fastq=rules.seqtk_fasta_to_fastq_full.output,
        index=rules.sencha_index.output
    output:
        coding_prot=os.path.join(translate_dir, "{samplename}_{rep}.{alphabet}_k{ksize}_ref{ref}_t{tablesize}_jacc{jaccard_thresh}.codingpep.fa"),
        coding_nucl=os.path.join(translate_dir, "{samplename}_{rep}.{alphabet}_k{ksize}_ref{ref}_t{tablesize}_jacc{jaccard_thresh}.codingnucl.fa"),
        noncoding_nucl=os.path.join(translate_dir, "{samplename}_{rep}.{alphabet}_k{ksize}_ref{ref}_t{tablesize}_jacc{jaccard_thresh}.noncoding.fa"),
        low_complexity_prot=os.path.join(translate_dir, "{samplename}_{rep}.{alphabet}_k{ksize}_ref{ref}_t{tablesize}_jacc{jaccard_thresh}.lowcomplexprot.fa"),
        low_complexity_nucl=os.path.join(translate_dir, "{samplename}_{rep}.{alphabet}_k{ksize}_ref{ref}_t{tablesize}_jacc{jaccard_thresh}.lowcomplexnucl.fa"),
        parquet=os.path.join(translate_dir, "{samplename}_{rep}.{alphabet}_k{ksize}_ref{ref}_t{tablesize}_jacc{jaccard_thresh}.parquet"),
        json=os.path.join(translate_dir, "{samplename}_{rep}.{alphabet}_k{ksize}_ref{ref}_t{tablesize}_jacc{jaccard_thresh}.json"),
    log: 
        os.path.join(logs_dir, "sencha_translate", "{samplename}_{rep}.{alphabet}_k{ksize}_ref{ref}_t{tablesize}_jacc{jaccard_thresh}.translate.log") #2>{log} err ("missing PEPTIDES file")
    benchmark: 
        os.path.join(logs_dir, "sencha_translate", "{samplename}_{rep}.{alphabet}_k{ksize}_ref{ref}_t{tablesize}_jacc{jaccard_thresh}.translate.benchmark")
    params:
        alphabet=lambda w: alpha_abbreviations[w.alphabet],
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: attempt *10000,
        runtime=60
    wildcard_constraints:
        ref="\w+",
        alphabet="\w+",
        ksize="\d+",
        jaccard_thresh="\d\.\d*",
    conda: os.path.join(envs_dir, "sencha-env.yml")
    shell:
        """
        sencha translate --verbose --peptides-are-bloom-filter --alphabet {params.alphabet} --peptide-ksize {wildcards.ksize} \
        --jaccard-threshold {wildcards.jaccard_thresh} --noncoding-nucleotide-fasta {output.noncoding_nucl} \
        --low-complexity-nucleotide-fasta {output.low_complexity_nucl} --coding-nucleotide-fasta {output.coding_nucl} \
        --low-complexity-peptide-fasta {output.low_complexity_prot} --parquet {output.parquet} --json-summary {output.json} \
        {input.index} {input.fastq} > {output.coding_prot} 2> {log}
        """

rule summarize_translate_results:
    input: 
        translate_json=translate_json_targets,
        species_metadata=lambda w: config["species_metadata"]
    output: os.path.join(translate_dir, "{run_name}_summary.csv")
    log: os.path.join(logs_dir, "aggregate_json", "{run_name}.aggregate_json.log")
    benchmark: os.path.join(logs_dir, "aggregate_json", "{run_name}.aggregate_json.benchmark")
    conda: os.path.join(envs_dir, "pandas-env.yml")
    shell:
        """
        python scripts/aggregate_sencha_json_summaries.py {input.translate_json} --output-csv {output} --species-metadata {input.species_metadata}
        """




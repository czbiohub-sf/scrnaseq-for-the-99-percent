REFSEQ_FASTA=/home/olga/data_lg/czbiohub-reference/ncbi/refseq/releases/refseq-release201--2020-07-21/nonredundant-protein/complete__nonredundant_protein.faa.gz
BUSCO_METAZOA=/home/olga/data_sm/immune-evolution/databases/busco/orthodb-v10/metazoa_odb10/metazoa_odb10.fasta
OUTDIR_BASE=/mnt/ibm_sm/olga/immune-evolution/pipeline-results/human/kmermaid
BAM_BASE=/home/olga/data_sm/tabula-sapiens
WORK_DIR=/mnt/data_lg/olga/nextflow-kmermaid-bladder/
REFSEQ_DIAMOND=/mnt/ibm_sm/olga/immune-evolution/pipeline-results/bat/nf-predictorthologs--busco_metazoa/diamond/complete__nonredundant_protein_db.dmnd
REF_PROTEOME=/home/olga/data_lg/czbiohub-reference/uniprot/releases/2019_11/manually_downloaded/uniprot-reviewed_yes+taxonomy_Mammalia_9MAMM_40674.fasta.gz

bladder:
	nextflow run \
		nf-core/kmermaid \
		-r  olgabot/sourmash-sig-merge \
		-latest \
		-profile docker \
		 --bam ${BAM_BASE}/'**.bam' \
		 --reference_proteome_fasta ${REF_PROTEOME} \
		 --sketch_scaled 10 \
		 -with-tower \
		 -work-dir ${WORK_DIR} \
		 --outdir ${OUTDIR_BASE}/$@ \
		 -resume \
		 --tenx_min_umi_per_cell 2500 \
		 --molecules dayhoff \
		 --track_abundance \
		 --ksizes 27,30,33,36,39,42,45,48,51 \
		 --skip_compare \
		--skip_compute

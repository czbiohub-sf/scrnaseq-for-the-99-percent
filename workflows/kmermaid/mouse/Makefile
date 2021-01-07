REFSEQ_FASTA=/data_lg/olga/czbiohub-reference/ncbi/refseq/releases/refseq-release201--2020-07-21/nonredundant-protein/complete__nonredundant_protein.faa.gz
OUTDIR_BASE=/data_sm/olga/immune-evolution/pipeline-results/mouse/kmermaid
TGZ_BASE=/data_lg/olga/tabula-muris-senis/10x/per-tissue-softlinks/Bladder
WORK_DIR=/data_lg/olga/nextflow-kmermaid-bladder/
REF_PROTEOME=/data_lg/olga/czbiohub-reference/uniprot/releases/2019_11/manually_downloaded/uniprot-reviewed_yes+taxonomy_Mammalia_9MAMM_40674.fasta.gz

bladder:
	nextflow run \
		nf-core/kmermaid \
		-r olgabot/sourmash-sig-merge \
		-latest \
		--tenx_tgz ${TGZ_BASE}/'**.tgz' \
		--reference_proteome_fasta ${REF_PROTEOME} \
		--sketch_scaled 10 \
		-work-dir ${WORK_DIR} \
		--outdir ${OUTDIR_BASE}/$@ \
		-resume \
		--tenx_min_umi_per_cell 2500 \
		--molecules dayhoff \
		--track_abundance \
		--ksizes 27,30,33,36,39,42,45,48,51 \
		--skip_compare

BAM_DIR=/home/olga/data_sm/tabula-microcebus/analyses/kmermaid/tenx-tgz--alldata/10x-bams/immune-softlinks
TRANSLATE_FASTA=/home/olga/data_lg/czbiohub-reference/uniprot/releases/2019_11/manually_downloaded/uniprot-reviewed_yes+taxonomy_Mammalia_9MAMM_40674.fasta.gz

antoine:
	nextflow run \
                nf-core/kmermaid \
		-r olgabot/grep-per-cell-fastqs \
		-latest \
		-profile docker \
		--bam ${BAM_DIR}'/*LUNG*.bam' \
		--skip_compute \
		--skip_compare \
		--molecules dayhoff \
		--ksizes 27,30,33,51,54,57,60 \
		--sketch_scaled 10 \
		--track_abundance \
		--reference_proteome_fasta ${TRANSLATE_FASTA} \
		--outdir ~/data_sm/tabula-microcebus/analyses/kmermaid/tenx-tgz--alldata-v2-from-bams \
		-work-dir ~/data_sm/nextflow-intermediates/ \
		 --tenx_min_umi_per_cell 5000 \
		-with-tower \
		-resume

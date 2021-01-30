REFERENCE_PROTEOME_FASTA=/home/olga/data_lg/czbiohub-reference/uniprot/taxonomy_Mammalia_9MAMM_40674_plus_HLrhiFer5_Bat1K.fasta.gz
OUTDIR_BASE=/mnt/ibm_sm/olga/immune-evolution/pipeline-results/bat/kmermaid/
COMMON_FLAGS=-latest \
		-resume \
		-with-tower \
		--reference_proteome_fasta ${REFERENCE_PROTEOME_FASTA} \
		--tenx_min_umi_per_cell 2500 \
		--skip_compute \
		--skip_compare \
		--skip_multiqc \
		--sketch_scaled 10 \


redo_singlecell_fastas:
	nextflow run nf-core/kmermaid \
		-r olgabot/translate-outputs \
		-work-dir /data_sm/olga/nextflow-kmermaid-work/ \
		${COMMON_FLAGS} \
		--bam '/data_sm/olga/batlas/Bat_Bams/*LU*.bam' \
		--outdir ${OUTDIR_BASE}/$@


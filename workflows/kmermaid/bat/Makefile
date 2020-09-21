REFERENCE_PROTEOME_FASTA=/home/olga/data_lg/czbiohub-reference/uniprot/taxonomy_Mammalia_9MAMM_40674_plus_HLrhiFer5_Bat1K.fasta.gz
BAM_DIR=/mnt/ibm_sm/home/phoenix/batlas/BatBams/immune-tissue-softlinks
OUTDIR_BASE=/mnt/ibm_sm/olga/immune-evolution/pipeline-results/bat/kmermaid
WORK_DIR=/mnt/ibm_sm/home/phoenix/nextflow-work

run-1:
	nextflow run nf-core/kmermaid \
		-r phoenix-bam-input-fix \
		-latest \
		-profile docker \
		-resume \
		-work-dir ${WORK_DIR} \
		--molecules protein,dayhoff \
		--sketch_scaled 10 \
		--translate_peptide_ksize 11 \
		--reference_proteome_fasta ${REFERENCE_PROTEOME_FASTA} \
		--tenx_min_umi_per_cell 1000 \
		--bam '/mnt/ibm_sm/home/phoenix/batlas/BatBams/immune-tissue-softlinks/*.bam' \
		--track_abundance \
		--outdir ${OUTDIR_BASE}/$@ \